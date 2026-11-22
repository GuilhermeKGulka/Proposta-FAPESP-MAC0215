# PLANO — Gencan com CG em GPU (hessiana explícita)

> Documento de migração de contexto. Resume o estado do projeto, decisões, código, pendências e próximos passos.
> Itens marcados **[NÃO VERIFICADO]** foram escritos mas ainda não compilados/rodados. Itens **[ABERTO]** dependem de informação ou decisão pendente.
> Última atualização: 22/11/2026.

---

## 1. Objetivo

Levar para a GPU o **gradiente conjugado (`cgm`) do Gencan**, mantendo o método original (laço externo, busca linear, extrapolação, projeção, faces livres/`ifr`, critérios de parada do `cgm`) e trocando apenas o produto Hessiana-vetor por SpMV na GPU (hipSPARSE, API genérica) mais kernels HIP para as reduções e atualizações.

Escopo atual: **somente hessiana explícita** (`truehp = .true.`). A versão com aproximação secante (`truehp = .false.`) fica fora por enquanto.

Estilo de implementação: o mesmo dos benchmarks de "API genérica" do repositório (`Generic-SpMV-API-Benchmarks/CUTEst-Newton-CG-test` é o exemplo mais próximo): Fortran + hipfort + bindings `bind(c)` próprios para a API genérica do hipSPARSE (`hipsparse_generic.f90`), CMake com presets `amd`/`nvidia`, kernels em `.cpp` compilados com `hipcc`.

Hardware de desenvolvimento: AMD RX 9070 XT (`gfx1201`, ROCm, gfortran 16). Também há backend NVIDIA via hipfort, mas os testes recentes foram só em AMD.

---

## 2. Estrutura do repositório (relevante)

| Diretório | Papel |
|---|---|
| `gencan-cg-standalone/` | Gencan **original (CPU)**. Referência. **Não alterar** `gencan.f90`, `cgm.f90`, `lssit.f90`. |
| `gencan-cg-standalone-gpu/` | Trabalho atual. Contém a reescrita (ver §4). A tentativa antiga de GPU foi descartada. |
| `Generic-SpMV-API-Benchmarks/` | Benchmarks SpMV/CG/Newton-CG com API genérica (base de estilo e de `hipsparse_generic.f90`). |
| `Initial-Benchmarks/` | Benchmarks anteriores (API legada `csrmv`), Julia x Fortran, resultados de análise. |
| `externals/` | CUTEst, SIFDecode, ARCHDefs (submódulos), `mastsif/` (não versionado). |
| `env_cutest.sh` | `source ../env_cutest.sh` define `CUTEST`, `SIFDECODE`, `MASTSIF`, `MYARCH` (obrigatório até para o CMake). |

Resultados prévios de benchmark (contexto): CG sobre SuiteSparse, GPU AMD, matrizes onde ambos convergem: speedup mediana ≈ 7.6x, máximo ≈ 96x. Ganhos aparecem com NNZ grande; com n pequeno a latência de sincronização domina.

---

## 3. Por que recomeçar (erros da tentativa antiga `gencan-gpu`)

Registrados para não repetir:

1. `d` nunca voltava ao host (todas as saídas eram `return`, o `hipMemcpy` ficava em código inalcançável).
2. Faltava `hipblasSetPointerMode(handle, 1)` (escalares em device).
3. Inversão de direção errada: com `ptr > 0` o original faz `p = -p`; o código fazia `p = -z`.
4. `abox` calculado com o `d` antigo (no original usa o `d` já atualizado).
5. `iter` incrementado tarde demais, saídas antecipadas com `iter` errado.
6. Faltavam critérios do `cgm`: ângulo (cginfo=3), modelo quadrático (4), `samep` (6), reinício sem precondicionador; `precond` ignorado.
7. `atomicMin` em `double` via `unsigned long long` quebra com valores negativos.
8. `copy_hessian_to_gpu` refazia CSR, descritor, `bufferSize` e buffer a **cada chamada** do CG.
9. Matriz criada `n×n` com índices renumerados para `1..nfr` (funcionava por acaso).
10. Possíveis causas do `HIPSPARSE_STATUS_INVALID_VALUE` (código 3): `nnz = 0`, descritores de vetor com tamanho ≠ da matriz, buffer não nulo quando `bufferSize = 0`, preprocess incoerente.

---

## 4. Código existente (nova arquitetura)

Todos em `gencan-cg-standalone-gpu/` (fontes em `sources/`, testes em `tests/`).

### 4.1 Arquivos

| Arquivo | Origem | Conteúdo / estado |
|---|---|---|
| `sources/lssit.f90`, `cgm.f90`, `gencan.f90` | cópia **intacta** do `gencan-cg-standalone` | Referência CPU. Biblioteca CMake `gencan_cpu`. |
| `sources/hipsparse_generic.f90` | reaproveitado dos benchmarks | Bindings `hipsparseCreateCsr`, `CreateDnVec`, `SpMV_bufferSize/preprocess/SpMV`; `hs_alloc_buffer` devolve ponteiro **nulo** se `bufferSize = 0`. |
| `sources/hess_csr.f90` | novo | COO (do `evalh`) → CSR `nfr×nfr` na GPU + SpMV. Ver §4.2. Versão com **scatter em kernel**, validada em 21/11. |
| `sources/cgm_kernels.cpp` | novo | Kernels do CG (reduções, atualizações, scatter). Versão da fase 5, validada em hardware AMD em 21/11. |
| `sources/cgm_gpu.f90` | novo | `cgm_gpu` fiel ao `cgm`, contexto `cgm_gpu_ctx`. Fase 5 (1 sincronização por iteração) validada em 21/11. |
| `sources/gencan_gpu.f90` | cópia do `gencan.f90` original com mudanças mínimas | Módulo `modgencan_gpu`. Versão corrigida entregue. Ver §4.4. |
| `sources/cutestma_gpu.f90` | cópia do `cutestma.f90` | Só muda `use modgencan_gpu, only: gencan, genunc`. |
| `tests/test_hess_csr.f90` | novo | Etapa 1: CSR/SpMV contra `preparehv`+`computehv` reais. |
| `tests/test_cgm_gpu.f90` | novo | Etapas 3–4: `cgm_gpu` contra `cgm` original em cenários sintéticos. Tem **9 cenários** (7–9 com precondicionador) e imprime OK/FALHOU por cenário. |
| `run-cutest-gpu.sh` | escrito na conversa | **Já tem o conteúdo novo** (conferido em 21/11): roda `cutestma_cpu` e `cutestma_gpu` por problema → `table_cpu.txt`/`table_gpu.txt`, retomável via `done_gencan.txt`. Não existe `run_cutest_gpu.sh`. |
| `run-tests.sh` | escrito na conversa | Compila e roda `test_hess_csr` e `test_cgm_gpu`. Funciona: o alvo `test_hess_csr` voltou ao `CMakeLists.txt` em 21/11 (tinha sido removido em bff320c). |
| `run-cutest.sh`, `compile.sh` | tentativa antiga | Obsoletos: usam alvos `cutestma`/`toyprob`, que não existem com esses nomes (desde a fase 6 há `toyprob_cpu`). `README.md` também está desatualizado. Saem no passo 6 da fase 6. |
| `CMakeLists.txt`, `CMakePresets.json` | fase 6, passo 2 (22/11) | GPU opcional (`GENCAN_GPU` AUTO/ON/OFF) e CUTEst opcional. Presets `amd`/`nvidia` (GPU obrigatória) e `cpu` (sem HIP). Ver §4.5. |

### 4.2 `hess_csr_mod` (Hessiana CSR por face)

- Entrada: COO em índices `1..n` como vem do `evalh` (triângulo, diagonal + fora da diagonal espelhada) e `frvar(n)`.
- Renumeração **idêntica à de `preparehv`**: i-ésimo `.true.` de `frvar` → i (consistente com `ifr = pack(1..n, frvar)`). **Não se chama `preparehv`** no caminho GPU: ele destrói o COO original.
- Para cada entrada COO `e`: `pos1(e)` = posição de `(r,c)` no CSR, `pos2(e)` = posição de `(c,r)` (0 se diagonal); ambos 0 se a entrada sai da face. Atualizar valores = scatter simples.
- **Padrão só é refeito se a face muda** (ou `hne` muda, ou `hess_csr_invalidate`). Então destrói/recria matriz, vetores, buffer e preprocess. Se só os valores mudam, reaproveita descritor e buffer (como no Newton-CG).
- Guardas: `nnz = 0` ou `nfr = 0` → `y = 0` via `hipMemset`, sem chamar hipSPARSE.
- SpMV usa `d_x = p` e `d_y = hp` do contexto do CG (ponteiros passados no `hess_csr_init`).
- Scatter (versão nova): copia `hne` doubles COO para `d_hcoo` e `launch_scatter` escreve em `d_val`; `pos1/pos2` sobem só quando o padrão muda.
- `hess_csr_invalidate(h)`: força refazer o padrão. **Necessário nos testes sintéticos** (cada cenário tem COO novo com a mesma face); desnecessário no Gencan real (estrutura do COO do CUTEst é fixa por problema).

### 4.3 `cgm_gpu` (fiel ao `cgm`)

Interface igual à do `cgm`, trocando `lssit` por `ctx` (e passando `lssit_obj` para ler os dados do precondicionador). Antes de chamar: `hess_csr_update(ctx%h, ...)` e, se `precond`, `comphpre(...)` na CPU.

Pontos de projeto (todos validados nas fases 3–4, exceto onde indicado):

- **Tudo sobre as variáveis livres**: vetores e CSR com tamanho `nfr`. `x`, `lbnd`, `ubnd`, `lind`, `uind` são compactados por `ifr` no host **uma vez por chamada** e enviados (`lind/uind` como flags inteiros, para `abox` exato).
- **Sem inverter `p` na GPU.** O original faz `p = -p` quando `ptr > 0`. `k_dots1` calcula `amax` para `+p` e `-p`; o host escolhe com `sgn = ±1` (`ptd`, `gtp`, `coef = sgn*alpha` carregam o sinal; `p·hp` e `p·p` independem do sinal; o `coef` do `p` seguinte usa `sgnprev`).
- **Reduções**: kernels próprios, redução por bloco + redução final no último bloco (contador + `__threadfence`), um lançamento, determinístico, sem `atomicMin/atomicAdd` em double. `BS = 256`, árvore em memória compartilhada (independe de wavefront 32/64), grade ≤ 1024 blocos com grid-stride. Contador `unsigned` reinicia sozinho (`atomicInc`).
- **Ordem dos testes idêntica ao `cgm`**: parada 0 e 8; `pthp`; ângulo (cginfo 3); `q`; atualização; `iter+1`; `alpha == amax` (1); `abox <= 0.1` (5); `samep` (6); progresso de `q` (4); cginfo 7.
- **Reinício sem precondicionador** (`restarted`, `itertmp`, `cycle outer`) já estruturado.
- **Precondicionador (m = p = 0)**: `z = r/D + c1·u + c2·s`, com `D = plspg + pdiag`, `u = psmdy`, `c1 = s·r/sty`, `c2 = u·r/sty − psmdyty·(s·r)/sty²`. `comphpre` continua na CPU (O(n), já validado). `ztr = Σr²/D + c1·(u·r) + c2·(s·r)` usa 3 escalares que `k_update` calcula **na mesma passada** (sem sincronização extra); `k_dots0` os calcula no início/reinício. `znorm2` não é necessário (usa-se `p·p`). Se `corr` é falso, `s` e `u` são zerados para evitar `0·NaN`.
- **Fase 5 (validada em 21/11)**: uma sincronização por iteração via **especulação**. Após decidir o passo de k, o host lança `update(k)` → `update_p_dev(k+1)` (calcula `beta` na GPU lendo `red`) → SpMV → `dots1(k+1)` e só então lê `red` (`d_red`, 16 doubles: `[0..5]` de `k_update`/`k_dots0`, `[8..14]` de `k_dots1`). Se o CG para, o trabalho especulativo (só escreve em `p` e `hp`) é descartado; `d` e `r` já estão corretos. Custo: no máximo uma iteração desperdiçada por chamada.
  - Risco conhecido: o `beta` no kernel pode usar FMA e diferir 1 ulp do host; pode mudar `iter` em 1 em casos sensíveis.

Layout de `d_red` (fase 5):

| Índice (0-based) | Origem | Conteúdo |
|---|---|---|
| 0 | `k_update` | `r·r` |
| 1 | `k_update` | nº de violações de `samep` |
| 2, 3, 4 | `k_update` / `k_dots0` | `Σr²/D`, `s·r`, `u·r` |
| 5 | `k_update` | `abox` |
| 8..14 | `k_dots1` | `p·hp`, `p·d`, `g·p`, `p·r`, `p·p`, `amax(+p)`, `amax(−p)` |

### 4.4 `gencan_gpu.f90` (mudanças em relação ao original)

- Módulo `modgencan_gpu`; `use hess_csr_mod`, `use cgm_gpu_mod` no lugar de `use modcgm, only: cgm`.
- `genunc` e `gencan`: `type(cgm_gpu_ctx) :: ctx` **local**; `lssitini`, `cgm_gpu_init(ctx, n, hnnzmax, .true.)`, chama `genuncb`/`gencanb` com `ctx` como **último argumento**, `cgm_gpu_free`, `lssitend`. `stop` se `.not. truehp`.
- `genuncb` e `gencanb`: novo argumento `type(cgm_gpu_ctx), intent(inout) :: ctx` e `logical :: changed, pc`.
- Em `genuncb` e em `newtonls`: **sem `preparehv`**; no lugar:
  `call hess_csr_update(ctx%h, lssit%hne, lssit%hrow, lssit%hcol, lssit%hval, frvar, changed)`;
  `pc = .false.; if (precond) pc = (inform .eq. 0)`;
  `call cgm_gpu(ctx, n, nfr|n, ifr, x, lind, lbnd, uind, ubnd, g(ifr), cgdel, cgeps, maxitcg, pc, d, cgiter, iprintcg, ierr, cginfo, lssit)`.
- O `comphpre` e o resto (busca linear, extrapolação, projeção) ficam na CPU, intactos.
- (O `pc` também corrige o original, que lê `inform` sem inicializar quando `precond = .false.`.)

### 4.5 CMake (alvos esperados)

Desde a fase 6, passo 2 (22/11/2026):
- **`GENCAN_GPU`** (`AUTO` por padrão):
  - `AUTO` procura o backend de `HIP_PLATFORM` (amd: hipfort, hip, hipsparse, hipblas, hipcc e compilador C++; nvidia: hipfort nvptx, CUDA Toolkit, hipsparse/hipblas para CUDA e hipcc). Se faltar algo, avisa o que falta e compila só a CPU.
  - `ON` aborta listando o que falta.
  - `OFF` nem procura.
  - A lógica de busca é a dos `Generic-SpMV-API-Benchmarks`, testada em ROCm e CUDA. Só deixou de ser fatal.
- **Presets:** `amd` e `nvidia` usam `GENCAN_GPU=ON`; `cpu` usa `OFF` e gera `build/cpu`.
- **CUTEst opcional:** sem as variáveis do `env_cutest.sh`, os alvos `cutestma_*` não são criados.
- **Módulos:** os de `gencan_cpu` ficam em `mod/cpu` e os da parte GPU em `mod/gpu`, preparando a segunda compilação dos mesmos fontes com `-DGENCAN_GPU` (D2(b)).
- **Pré-processamento:** `lssit/cgm/gencan` são compilados com `Fortran_PREPROCESS ON` (`-cpp`).
- **Verificado em 22/11:**
  - o build `cpu` não liga nada de HIP (`ldd`);
  - a saída de `toyprob_cpu` e `cutestma_cpu` é **idêntica** à da referência compilada pelos scripts dela, nos 7 problemas da `pequena.txt` (fora a linha de tempo);
  - nvidia AUTO sem CUDA cai para só CPU, e nvidia ON aborta;
  - amd AUTO com `hip` ausente (simulado) cai para só CPU;
  - preset `amd`: testes 5/5 e 9/9 e ROSENBR igual em CPU e GPU.
  - O build NVIDIA com CUDA ainda não foi testado com o código novo **[NÃO VERIFICADO]**.

```
toyprob_cpu  : sources/toyprob.f90                            (liga gencan_cpu; sempre)
gencan_cpu   : sources/lssit.f90 cgm.f90 gencan.f90          (referência, inalterados)
hsgen        : sources/hipsparse_generic.f90
gpucore      : sources/hess_csr.f90 cgm_gpu.f90              (liga hsgen, hipfort_backend, gencan_cpu)
kernels      : sources/cgm_kernels.cpp  (hipcc, -DKERNEL_FLAGS="--offload-arch=gfx1201")
gencan_gpu   : sources/gencan_gpu.f90                        (liga gpucore)
test_hess_csr, test_cgm_gpu : tests/*.f90  (ligam gpucore gencan_cpu kernels hsgen hipfort_backend)
cutestma_cpu : sources/cutestma.f90    + ${DECODED_SRCS} + ${CUTEST_LIB}   (liga gencan_cpu)
cutestma_gpu : sources/cutestma_gpu.f90 + ${DECODED_SRCS} + ${CUTEST_LIB}  (liga gencan_gpu gpucore kernels hsgen hipfort_backend)
```

**Armadilha**: `DECODED_SRCS` e `CUTEST_LIB` só existem depois do bloco CUTEst (sifdecoder) do `CMakeLists.txt`. Os `add_executable(cutestma_*)` precisam vir **depois** dele, senão o link falha com `undefined reference to elfun_/group_/range_`.

---

## 5. Decisões tomadas (resumo)

1. Reescrever do zero a partir do Gencan/`cgm` **originais da CPU**; a versão GPU antiga não é base.
2. Hessiana **CSR `nfr×nfr`** (submatriz das livres), renumeração igual a `preparehv`, sem chamar `preparehv`.
3. Padrão CSR/descritor/buffer/preprocess recriados **só quando a face muda**; valores por scatter.
4. `bufferSize = 0` ⇒ buffer **nulo**; guarda para `nnz = 0` / `nfr = 0`.
5. Kernels de redução próprios (sem hipBLAS no `cgm`); `hipblasSetPointerMode` deixa de ser necessário.
6. Sinal de `p` tratado por `sgn` (sem kernel de inversão).
7. Compactação de `x/lbnd/ubnd/lind/uind` por `ifr` uma vez por chamada.
8. Precondicionador: `comphpre` na CPU; aplicação fundida nas reduções existentes.
9. Fase 5: 1 sincronização por iteração via especulação; scatter COO→CSR em kernel.
10. Validar sempre contra o `cgm`/`gencan` **originais** (mesmo `iter`, `cginfo`, `fcnt/gcnt/hcnt`, `f`, `gpsupn`).
11. Fase 6 (22/11/2026; detalhes no `docs/FASE6.md` §4.0):
    - código único com `#ifdef GENCAN_GPU`;
    - por enquanto, CPU e GPU escolhidas na compilação (duas bibliotecas a partir dos mesmos fontes);
    - estado da GPU no `lssit_type`;
    - CSR montado por espelhamento + `coo2csr` do orientador dentro do `preparehv`, a medir contra o scatter;
    - trabalho no branch `fase6`.
12. Fase 6, D4 medido (22/11/2026):
    - No modo GPU, o `preparehv` **não faz a compactação da CPU**. Ele parte do COO original, monta o padrão (espelho + `coo2csr` com a origem de cada entrada → `pos1/pos2`) só quando a estrutura ou a face mudam, e a cada passo envia os valores e espalha por kernel.
    - Resultado: 1,7 ms por passo, contra 55 ms da versão com a compactação e 1,05 ms da fase 5.
    - Toda a memória da GPU é alocada uma vez.
    - Em aberto: `coo2csr` in-place × contagem na troca de face (46 contra 14 ms), e comparação completa da estrutura (+0,65 ms/passo). Ver `docs/FASE6.md` §4.0.
    - Decidido em 22/11: contagem fora do lugar, como os dois primeiros laços do `coo2csr`, mas sem a permutação por ciclos. A face nova custa 16 ms (fase 5: 14). A versão com o `coo2csr` está no branch `fase6-coo2csr`.
    - A comparação completa da estrutura fica, por segurança (~0,6 ms por passo, limitada pela banda de memória). Será reavaliada com o perfil do Gencan real.

---

## 6. Estado dos testes

| Etapa | Teste | Resultado |
|---|---|---|
| 1. COO→CSR + SpMV | `test_hess_csr` | OK. Na fase 5 (21/11): erro relativo 0 a 1.5e-16 nos 5 trials. Fase 4: 9e-17 a 2.5e-16. Desde 21/11 o teste roda sem `hess_csr_invalidate` (o padrão COO é fixo, como no CUTEst) e confere o `changed` esperado: T, T, **F**, T, T. O trial 3 cobre "mesmo padrão, valores novos" (só o scatter): 5/5 OK. |
| 3. `cgm_gpu` sem precond | `test_cgm_gpu` cenários 1–6 | OK: `iter` e `cginfo` idênticos à CPU. Erro de `d` de 2e-16 a 1.2e-15 na fase 4 e de 0 a 1.5e-15 na fase 5 (21/11). Requer `call hess_csr_invalidate(ctx%h)` antes de cada `hess_csr_update` no teste (COO novo a cada cenário). |
| 4. precondicionador | `test_cgm_gpu` cenários 7–9 | OK (21/11): 7 = SPD com correção BFGS (`sameface = .true.`), 8 = bounds próximos com precond, 9 = SPD sem correção (`sameface = .false.`). `iter` e `cginfo` idênticos, erro de `d` ≤ 2.5e-16. Critério do teste: mesmo `iter`, mesmo `cginfo`, erro relativo de `d` ≤ 1e-12. Total: 9/9 OK. Obs.: a semente aleatória não é fixada, então `nfr` e `iter` variam entre execuções. |
| Integração Gencan | `cutestma_cpu` vs `cutestma_gpu`, 7 problemas | `iter`, `fcnt`, `hcnt` idênticos; `gcnt` da GPU = CPU + 1; `f` difere no último dígito; **`BQPGAUSS`: 97 (CPU) vs 118 (GPU) iterações**. |
| 5. fase 5 | todos os anteriores | **Validada em 21/11/2026** (detalhes abaixo). Duas correções foram necessárias só para compilar: o atributo `target` num componente de tipo em `hess_csr.f90`, e o argumento `ls` que faltava na chamada de `cgm_gpu` em `test_cgm_gpu.f90`. Resultados em `gencan-cg-standalone-gpu/resultados/table_{cpu,gpu}_fase5_{nopre,pre}.txt`. |

Resultados CUTEst (CPU vs GPU, fase 4):

| Problema | n | iter CPU/GPU | gcnt CPU/GPU | CPU(s) CPU/GPU |
|---|---|---|---|---|
| ROSENBR | 2 | 23/23 | 24/25 | 0.0003 / 0.157 |
| BDQRTIC | 5000 | 10/10 | 11/12 | 0.016 / 0.140 |
| CRAGGLVY | 5000 | 14/14 | 15/16 | 0.023 / 0.152 |
| EDENSCH | 2000 | 12/12 | 13/14 | 0.005 / 0.127 |
| BQPGAUSS | 2003 | **97/118** | 98/120 | 0.576 / 2.415 |
| TORSION1 | 5476 | 46/46 | 47/48 | 0.071 / 0.266 |
| OBSTCLAE | 10000 | 38/38 | 39/40 | 0.129 / 0.331 |

**Interpretação do `gcnt` (+1)**: é a chamada extra a `evalg` do `comphpre` na inicialização. CPU tem `gcnt = iter + 1` (⇒ `precond = .false.`), GPU `gcnt = iter + 2` (⇒ `precond = .true.`). Os dois `cutestma` estão com `precond` **diferente**. Confirmado em 21/11: `cutestma.f90:80` tem `.false.` e `cutestma_gpu.f90:80` tem `.true.`. A divergência do `BQPGAUSS` provavelmente é o mesmo efeito (hipótese, **[ABERTO]**).

**Fase 5, mesmo `precond` nos dois lados (21/11/2026).** Os drivers agora usam `precond = .true.`, igual ao `cutestma.f90` de referência.

| Problema | iter CPU/GPU, precond F | iter CPU/GPU, precond T | `fcnt/gcnt/hcnt` |
|---|---|---|---|
| ROSENBR | 23/23 | 23/23 | idênticos |
| BDQRTIC | 10/10 | 10/10 | idênticos |
| CRAGGLVY | 14/14 | 14/14 | idênticos |
| EDENSCH | 12/12 | 12/12 | idênticos |
| TORSION1 | 46/46 | 46/46 | idênticos |
| OBSTCLAE | 38/38 | 38/38 | idênticos |
| BQPGAUSS | **97/90** | **89/118** | diferem |

- Com o mesmo `precond`, o `gcnt` bate em todos os problemas. A hipótese do "+1" está confirmada.
- Fora do BQPGAUSS, `f` concorda em pelo menos 13 dígitos (critério 1e-6 da proposta atendido com folga), e `iter/fcnt/gcnt/hcnt` são idênticos.
- A GPU da fase 5 reproduz a da fase 4 iteração por iteração (BQPGAUSS = 118 nas duas). A especulação não mudou a trajetória, e o risco de FMA no `beta` não apareceu.

**BQPGAUSS: arredondamento, não bug.** Comparei a trajetória com `iprint = 10`, numa cópia temporária fora do repositório. A CPU e a GPU são idênticas até a iteração externa 10, CG iteração 16, onde o `rnorm` passa a diferir (1.7911 contra 1.7893). O `rnorm` oscila muito nesse trecho (1.8 → 7.4 → 2.8), sinal de perda de ortogonalidade num sistema mal condicionado.

Controle: rodei o `cgm` da **CPU** com os produtos internos somados em ordem reversa. A divergência aparece **no mesmo ponto** (1.7911 contra 1.7903), e o número final de iterações muda também na CPU: 97 sem precondicionador e 87 com (contra 89 da ordem normal). Ou seja, a mudança na ordem das somas que a redução em árvore da GPU introduz basta para separar as trajetórias neste problema. O `f` final concorda até 1e-12 em todos os casos, e todos convergem (SC 0).

**Tempos**: a GPU é mais lenta nesses n pequenos. Parte (~0.12 s fixos) é custo de `cgm_gpu_init` (handles do hipSPARSE, alocações) dentro da região medida. O `cutestma` mede com `cpu_time`, que não é tempo de parede com GPU.

**Fase 6, até o passo 4 (22/11/2026)** — histórico completo no `docs/CHANGELOG_GPU.md`:
- **Testes:**
  - `test_hess_csr`: 6/6. É o `preparehv` com GPU contra o `preparehv` + `computehv` da CPU, mais a medição D4.
  - `test_cgm_gpu`: 9/9. É o mesmo `cgm` com `usegpu = .false.` e com a GPU.
  - Os testes não precisam mais de `invalidate`: o `lssit_gpu` detecta sozinho a mudança de estrutura.
- **Build `cpu`:** `compare-cpu-ref.sh` dá 0 diferenças contra a referência.
- **CUTEst `pequena.txt`, `precond = .true.`:** CPU e GPU com **todos os campos idênticos aos da fase 5**, fora o tempo, inclusive o BQPGAUSS (118). Tabelas em `resultados/table_{cpu,gpu}_fase6_pre.txt`.
- **Padrões novos por chamada de `preparehv`:** 1 nos problemas sem caixa. Nos com caixa, quase todo passo: BQPGAUSS 109/111, TORSION1 23/23, OBSTCLAE 21/21.

---

## 7. Pendências e próximos passos (ordem sugerida)

### 7.1 Imediato (validação das fases 4–5)

Situação em 21/11/2026: **todos os itens do §7.1 estão feitos**, e a fase 5 está validada (§6). Em 22/11/2026 o plano da fase 6 foi aprovado (§7.5, `docs/FASE6.md`), e os passos 1–4 foram feitos no mesmo dia, no branch `fase6`.

Próximos passos da fase 6:
- passo 5: remover `gencan_gpu.f90`/`cutestma_gpu.f90`, que dependem de confirmação;
- (c): medir em problemas CUTEst grandes, com uma lista filtrada por `n` mínimo. **Avisar antes, porque é experimento longo**;
- passo 6: limpeza;
- passo 7: validação final, incluindo NVIDIA e a variante `precond = .false.`.

1. **[FEITO]** **Igualar `precond`** nos dois `cutestma` e rodar duas rodadas: `false/false` e `true/true` (`grep -n "precond *=" sources/cutestma*.f90`). Conferir `gcnt` iguais e `BQPGAUSS`. Se ainda divergir com `precond` igual: rodar com `iprint = 3` dos dois lados e comparar a trajetória.
2. **[FEITO, com 9 cenários]** **Completar `tests/test_cgm_gpu.f90` com 8 cenários**: 7 = SPD com precond, 8 = SPD com bounds e precond. Para cada cenário com precond: sortear `s`, calcular `y = H s` (H simétrica a partir do COO), `sty`, `sts`, chamar `comphpre(n, s, y, sty, sts, 1d-16, 1d16, .true., inform, ls)` antes do `cgm`/`cgm_gpu`; passar `prec(s)` no lugar de `.false.` e `ls` no final da chamada `cgm_gpu`. `nscen = 8`, novo array `prec`.
3. **[FEITO]** **Validar a fase 5** (arquivos já aplicados em 3e3f3d6, base da fase 4 já guardada em `resultados/table_*_fase4.txt`). Antes, devolver o alvo `test_hess_csr` ao `CMakeLists.txt` para que `run-tests.sh` volte a funcionar. Rodar `run-tests.sh` e o CUTEst; exigir os mesmos `iter`/`cginfo`. Diferença de 1 em `iter` ⇒ suspeitar de FMA no `beta` do kernel (decidir: aceitar ou calcular `ztr`/`beta` de forma idêntica ao host).
4. **[FEITO]** Corrigir o teste `test_hess_csr` para cobrir também "mesmo padrão, valores novos" (trial 3 sem `invalidate`).

### 7.2 Medição justa de desempenho

5. Medir com **tempo de parede** (`system_clock`) em `cutestma_cpu` e `cutestma_gpu`.
   - **Limite de tempo por execução:** o `run-cutest-gpu.sh` usa `ulimit -St 600`, que é tempo de CPU. Nos benchmarks de desempenho ele deve ser retirado ou bem maior, porque a vantagem da GPU cresce justamente nos problemas que demoram (observação do usuário, 22/11/2026). A medição (c) da fase 6 ainda rodou com 600 s.
6. Tirar `cgm_gpu_init` da região medida (contexto criado uma vez, ou medir separadamente), ou reportar "tempo total" e "tempo do laço" separados.
7. Benchmark com n grande: `test_cgm_gpu` com `n = 300000`, `noff = 900000` (antes/depois da fase 5); problemas CUTEst grandes (ex.: `YATP1CLS`, `YATP2CLS`, `FMINSURF`, `POWER`, `CURLY*`, `NCB20*`, séries `TORSION*`, `OBSTCL*`, `JNLBRNG*`, `DEGTRID*`). Lista completa em `problems.txt`.
8. Perfilar (rocprof): quanto vai em H2D por iteração de Newton (COO + compactações), `hess_csr_update`, SpMV, kernels e sincronizações.

### 7.3 Otimizações futuras (só depois de medir)

- Manter `x/lbnd/ubnd/lind/uind` compactados na GPU entre chamadas quando a face não muda.
- Reduzir cópias H2D por iteração de Newton (hoje: COO `hne` doubles + compactações + `g`).
- Fundir kernels (`update_p` + `dots`), ou usar HIP Graphs para as 4 chamadas por iteração.
- Reaproveitar o contexto entre chamadas de `gencan`/`genunc` quando possível.
- Precondicionador completo na GPU só se `comphpre` aparecer no perfil.
- Avaliar `hipsparseSpMV` com algoritmos alternativos (`HIPSPARSE_MV_ALG_*`) e preprocess (`nopre`).
- Verificar o backend NVIDIA (`cmake --preset nvidia`) com a API genérica.

### 7.4 Extensões de escopo (a decidir)

- `truehp = .false.` (secante): fora de escopo agora.
- Restrições/`m + p ≠ 0` (Algencan): o precondicionador GPU hoje exige `m = p = 0`; `cgm_gpu` faz `stop` caso contrário.
- Caso `nfr = 0` em `cgm_gpu`: o `hess_csr` trata o vazio, mas conferir que o Gencan nunca chama o `cgm` com `nfr = 0` (no original só chama `newtonls` com variáveis livres).

### 7.5 Pedidos externos e documentação

**Ordem combinada em 21/11/2026:** primeiro terminar o §7.1, até a fase 5 validada, na estrutura atual. Só depois começa a fase 6, descrita abaixo, que deve resultar na implementação provavelmente final.

#### Fase 6: código único CPU/GPU (pedido do orientador)

O pedido:
- O Gencan GPU deve **coexistir com o Gencan CPU no mesmo código**, sem um conjunto paralelo de arquivos novos.
- O mesmo pacote serve para quem quer GPU e para quem não quer. Sem HIP disponível, a compilação **não pode falhar**; a versão gerada é só CPU.
- **A estrutura do código GPU deve acompanhar a do Gencan original e funcionar, ao máximo, por analogia com ele.** Isso vale para os mesmos módulos, as mesmas rotinas e o mesmo fluxo, e não só para guardar a Hessiana no `lssit`. Adaptações que melhoram o desempenho são aceitas, mas cada uma precisa ser justificada e entrar no changelog.
- Reaproveitar o código e a estrutura que já existem:
  - um único `cutestma` (sem `cutestma_gpu.f90`);
  - o `lssit` continua guardando a Hessiana;
  - pouco código novo.
- A conversão COO → CSR deve ser uma adaptação da rotina `coo2csr` do orientador (in-place, O(nnz), sabidamente eficiente). A cópia está em `docs/coo2csr_orientador.f90`.

**Plano: `docs/FASE6.md`** (proposto em 21/11/2026 e aprovado em 22/11/2026). O ponto principal é que o "conflito" do `coo2csr` não existe no fluxo do Gencan: `evalh` reescreve o COO a cada passo de Newton, logo antes do `preparehv`. O estado da GPU passa a morar no `lssit`. Decisões (detalhes no FASE6 §4.0):
- D1: pré-processador, `GENCAN_GPU` AUTO/ON/OFF.
- D2: **provisoriamente** só em tempo de compilação. Os mesmos fontes são compilados em duas bibliotecas (CPU e GPU) e geram `cutestma_cpu`/`cutestma_gpu`. O argumento `usegpu` em `genunc`/`gencan` depende de o orientador aceitar a mudança de assinatura. **[ABERTO: perguntar ao orientador]**
- D3: estado da GPU no `lssit_type`.
- D4: espelhar + `coo2csr` in-place + enviar o CSR, **medindo** contra o scatter da fase 5 e procurando melhorias.
- D5: código em `gencan-cg-standalone-gpu/`.
- D6: branch `fase6`.

**Para conversar com o orientador [ABERTO]:**
1. **Mudar a assinatura de `genunc`/`gencan`.**
   - Para escolher CPU/GPU em tempo de execução: argumento `usegpu`, que é o D2(a).
   - Para o usuário avisar que a estrutura da Hessiana devolvida pelo `evalh` é fixa: evitaria comparar `hrow/hcol` a cada passo, ~0,6 ms em `nnz = 1,4 M` (`docs/FASE6.md` §4.0).
   - As duas mudanças são do mesmo tipo e podem ser decididas juntas.
2. **Conversão COO → CSR** (explicação completa, com corretude e eficiência, em `docs/CONVERSAO_CSR.md`). Usamos contagem por linha com colocação fora do lugar, que mantém os dois primeiros laços do `coo2csr` dele. A permutação in-place por ciclos ficou cerca de 3× mais lenta na troca de face (46 contra 16 ms), e nos problemas com caixa a face muda em quase todo passo. A versão com o `coo2csr` está no branch `fase6-coo2csr`.

Perguntas originais (respondidas no plano acima):
- Como escolher CPU ou GPU: opção do CMake (por exemplo `-DUSE_GPU=ON`, ligada só se o HIP for encontrado) com pré-processador (`#ifdef` em `.F90`), ou módulos separados escolhidos na hora de linkar. E se a escolha também existe em tempo de execução (um parâmetro como `truehp`/`precond`).
- Onde fica o código único: `gencan-cg-standalone/` continua sendo a referência intocada que serve de linha de base. O código unificado provavelmente nasce em `gencan-cg-standalone-gpu/`.
- Como adaptar o `coo2csr` (combinado em 21/11: este conflito será resolvido na fase 6). A rotina é in-place, mas o caminho GPU atual depende de manter o COO original: `preparehv` não é chamado justamente porque destrói o COO. Além disso, ela não filtra a face, não renumera e não espelha o triângulo, e a ordem das colunas dentro de cada linha não é garantida. É preciso decidir se a conversão roda sobre uma cópia, se reaproveita a renumeração do `preparehv`, e quanto do `hess_csr` (posições `pos1/pos2`, padrão refeito só quando a face muda) sobrevive.
- Critério de aceitação: o mesmo das fases 4–5. Os resultados da GPU unificada devem ser idênticos aos da fase 5 validada, que é a linha de base para a refatoração.

#### Documentação e changelog (redução da dívida cognitiva)

- `docs/CHANGELOG_GPU.md` (a criar **na fase 6**): registra, mudança por mudança, a diferença entre o Gencan de referência (`gencan-cg-standalone/`) e o Gencan GPU. Para cada mudança: arquivo/rotina, o que mudou, por quê e como foi validada.
  - Ele é escrito junto com a refatoração da fase 6, e não antes, porque a estrutura atual vai mudar.
  - Na fase 6, cada commit que mexe no código GPU atualiza o changelog.
- Documento de arquitetura (`docs/GENCAN_GPU.md`, depois da fase 6): quem roda onde, os kernels, o que fica idêntico ao original, as decisões (§5), os erros encontrados (§3 e §8), o protocolo de teste (§6, §9) e a tabela de resultados CPU × GPU.

---

## 8. Armadilhas já encontradas (para não repetir)

1. **Fortran não diferencia maiúsculas**: `ctx%d_D` e `ctx%d_d` são o mesmo nome (e dummies `D`/`d` na mesma interface). Resolvido renomeando a diagonal do precondicionador para `d_dg` (`dgf`, `dg`). Há um script de verificação rápida: listar identificadores únicos ignorando caixa e procurar colisões.
2. **Ordem do CMake**: `DECODED_SRCS`/`CUTEST_LIB` precisam estar definidos antes dos executáveis que os usam.
3. **`hess_csr_update` só refaz o padrão se a face/`hne` mudam**: testes sintéticos com COO novo precisam de `hess_csr_invalidate`. Sem isso, o cenário 6 do `test_cgm_gpu` dava erro 1.0.
4. **Edição manual do `gencan.f90`**: houve reticências literais (`call gencanb(..., evalstp,...)`), `ctx` não declarado e bloco duplicado em `genunc`. O `ctx` é local em `genunc`/`gencan`; só `genuncb`/`gencanb` o recebem como argumento.
5. **`gcnt` não bate se `precond` for diferente** entre os `cutestma`.
6. Scripts: o `run-cutest-gpu.sh` atual já roda CPU e GPU. `run-cutest.sh` e `compile.sh` são obsoletos, e `run-tests.sh` depende do alvo `test_hess_csr`, que sumiu do CMake.
7. O script marca "falhou compilar" sem mostrar o motivo: ver `build/<preset>/logs/<PROBLEMA>.log`.
8. **gfortran não aceita `target` em componente de tipo derivado.** O atributo vem do dummy: um argumento `type(...), target :: h` torna `c_loc(h%comp)` válido.
9. `run-cutest-gpu.sh` pula tudo o que está em `done_gencan.txt`. Para rodar de novo a mesma lista, esvazie esse arquivo e as tabelas (deixe só o cabeçalho).
10. **Problemas sensíveis a arredondamento** (ex.: BQPGAUSS): só mudar a ordem das somas já altera a trajetória. Antes de tratar uma divergência como bug, faça o controle CPU contra CPU com a ordem das somas invertida (§6).

---

## 9. Como rodar

```bash
cd gencan-cg-standalone-gpu
source ../env_cutest.sh

# Testes sintéticos (etapas 1, 3, 4)
cmake --preset amd -DKERNEL_FLAGS="--offload-arch=gfx1201"
cmake --build --preset amd --target test_hess_csr test_cgm_gpu
./build/amd/test_hess_csr      # esperado: erro rel ~1e-16; trial 5 com nfr=0
./build/amd/test_cgm_gpu       # esperado: mesmo iter/cginfo que a CPU; erro de d 1e-16..1e-14

# Problemas CUTEst, CPU x GPU
printf "ROSENBR\nBDQRTIC\nCRAGGLVY\nEDENSCH\nBQPGAUSS\nTORSION1\nOBSTCLAE\n" > pequena.txt
PRESET=amd ./run-cutest-gpu.sh pequena.txt      # roda CPU e GPU; pula o que esta em done_gencan.txt
paste -d'\n' table_cpu.txt table_gpu.txt | less -S

# Um problema na mao
cmake --preset amd -DPROBLEM=ROSENBR
cmake --build --preset amd --target cutestma_cpu cutestma_gpu
(cd build/amd && ./cutestma_cpu | tail -15 && ./cutestma_gpu | tail -15)

# So CPU, sem HIP/ROCm/CUDA (e sem CUTEst, se o env nao estiver carregado)
cmake --preset cpu && cmake --build --preset cpu && ./build/cpu/toyprob_cpu
```

Critério de aceitação geral: `iter`, `fcnt`, `gcnt`, `hcnt`, `f`, `gpsupn` iguais (ou muito próximos) entre CPU e GPU **com o mesmo `precond`**. Diferenças de arredondamento nas somas podem separar trajetórias em problemas mal condicionados; se acontecer em um bem condicionado, tratar como bug.

---

## 10. Convenções para continuar o trabalho (Claude Code)

As regras de trabalho completas estão no `CLAUDE.md`, na raiz do repositório: propor o plano antes de editar, uma tarefa por vez, commits pequenos, avisar antes de experimentos longos, não mudar escopo nem cronograma sem consultar. Ao fim de cada sessão, atualizar neste arquivo a data, o §6, o §7 e o §8.

- Não editar `gencan-cg-standalone/` (referência). A versão GPU deve diferir do original **só** nos pontos do §4.4.
- Todo novo kernel/rotina precisa de teste contra o código CPU original antes de entrar no Gencan.
- Variáveis Fortran: evitar nomes que só diferem por caixa.
- Reduções devem ser independentes de tamanho de wavefront (hoje `BS = 256`, árvore em shared memory).
- Mudanças de desempenho: medir antes e depois (tempo de parede, mesmo `precond` nos dois lados) e registrar nas tabelas deste arquivo.

---

## 11. Cronograma e critérios da proposta (`docs/Proposta.pdf`)

Cronograma de 12 meses (§5 da proposta). O projeto começou em **01/08/2026**, então o mês 1 é agosto de 2026 e o fim é julho de 2027. Em 21/11/2026 estamos no **mês 4**.

| Mês | Atividade | Situação (21/11/2026) |
|---|---|---|
| 1–2 (ago–set/26) | Estudo do Algencan (Newton truncado, CG), ambientação HIP, CMake/Git | Feito (benchmarks em `Initial-Benchmarks/` e `Generic-SpMV-API-Benchmarks/`) |
| 3–4 (out–nov/26) | CG com Hessiana explícita; testes de correção pequenos; integração ao Algencan | Em andamento, no prazo: integrado ao Gencan isolado, validando a fase 5 |
| 5 (dez/26) | Fechar a versão com Hessiana explícita; primeiros experimentos de desempenho CPU × GPU | Não iniciado. Previsto: fase 6 (código único) e medição justa (§7.2) |
| 6–7 (jan–fev/27) | Variante por quocientes incrementais (`truehp = .false.`) em GPU; testes de correção | Não iniciado |
| 8–10 (mar–mai/27) | Experimentos completos: porte, esparsidade, AMD × NVIDIA; (opcional) método direto esparso; profiling | Não iniciado |
| 11–12 (jun–jul/27) | Análise final, relatório técnico, código público no GitHub com instruções | Não iniciado |

A fase 6 (refatoração para código único) não está no cronograma da proposta. Pelo calendário, ela cabe no fim do mês 4 e no mês 5, junto com a medição justa. Também falta definir se e quando a integração passa do Gencan isolado para o Algencan completo (restrições, `m + p ≠ 0`).

Critérios de validação (§6.3 da proposta):
1. `f`, gradiente e multiplicadores de Lagrange iguais aos da CPU dentro de 1e-6 (ou do critério da literatura do problema).
2. Número de iterações externas idêntico entre CPU e GPU.
3. Speedup > 1 para problemas com pelo menos 1e5 variáveis.

Métricas exigidas (§6.1): tempo total até a convergência; tempo médio por iteração do CG, incluindo transferências; speedup `T_CPU/T_GPU`; memória da GPU (pico e média). Cada experimento é repetido, e reportam-se média e desvio padrão. Os problemas vão de 1e3 a 1e6 variáveis. Hardware: Ryzen 9 7900 + RX 9070 XT, e Ryzen 7 9700X + RTX 5070 Ti.
