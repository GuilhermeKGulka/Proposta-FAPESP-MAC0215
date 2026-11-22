# Fase 6 — Gencan CPU/GPU em código único (plano)

> Status: **plano aprovado em 22/11/2026**, com as decisões do §4.0. Proposta escrita em 21/11/2026.
> Base: pedido do orientador (PLANO §7.5) e o código da fase 5 validada (commit ddf53fd).

## 1. Objetivo

Um só Gencan que compila e roda com ou sem GPU e que, lido ao lado do original, funciona por analogia com ele: mesmos módulos, mesmas rotinas e mesmo fluxo. A GPU entra em poucos pontos bem marcados, e cada adaptação de desempenho é justificada no `docs/CHANGELOG_GPU.md`.

## 2. O que o código atual já mostra

1. **A diferença entre `gencan.f90` e `gencan_gpu.f90` já é pequena.** O `diff` atual tem quatro pontos:
   - criar e destruir o contexto (`cgm_gpu_init`/`cgm_gpu_free`) em `genunc`/`gencan`;
   - passar `ctx` para `genuncb`/`gencanb`;
   - trocar `preparehv` por `hess_csr_update`;
   - trocar `cgm` por `cgm_gpu`.

   O resto do arquivo é cópia idêntica, e por isso pode voltar a ser um arquivo só.
2. **O "conflito" do `coo2csr` in-place não existe no fluxo do Gencan.** A cada passo de Newton, o Gencan chama `evalh` (que reescreve `hrow/hcol/hval` no `lssit`) e logo depois `preparehv` (que já renumera e compacta o COO para a face, destruindo o original). O COO original nunca é reutilizado. O `hess_csr` atual guarda o COO só para uma otimização: refazer o padrão apenas quando a face muda e, fora isso, só espalhar os valores. Ela pode ser mantida de outra forma (§4.3).
3. **O ponto natural de analogia é o `lssit`.** No original:
   - `preparehv` prepara a Hessiana da face;
   - `computehv` faz `H v`;
   - `cgm` nunca vê a matriz, só chama `computehv`/`applyhpre`.

   Na GPU o papel é o mesmo: `preparehv` prepara o CSR na GPU, e o CG na GPU faz o SpMV.
4. **Código morto:** `cutest_iface.f90` e `gencan_gpu_context.f90` não são usados por nada. `compile.sh` e `run-cutest.sh` apontam para alvos que não existem mais.

## 3. Arquitetura proposta

| Arquivo | Hoje | Na fase 6 |
|---|---|---|
| `gencan.f90` | referência intacta + cópia `gencan_gpu.f90` | **Um só.** Muda só o necessário para ligar a GPU (ver D2). Não há `ctx`, porque o estado da GPU mora no `lssit`. |
| `lssit.f90` | referência intacta | Ganha o estado da GPU no `lssit_type` (`usegpu` + componente `gpu`). `lssitini`/`lssitend` alocam e liberam a GPU. `preparehv` faz, depois do que já faz hoje, o espelhamento + `coo2csr` + envio à GPU. |
| `lssit_gpu.f90` (novo) | `hess_csr.f90` | A Hessiana no device: tipo com os ponteiros, alocação, `coo2csr` adaptado, envio do CSR, descritores do hipSPARSE e o SpMV. Substitui o `hess_csr.f90`. |
| `cgm.f90` | referência intacta | Um desvio no início: `if (lssit%usegpu) then; call cgm_gpu(...); return; end if`. O resto fica intacto. |
| `cgm_gpu.f90` | idem | Mantido, porque a organização com especulação e uma sincronização por iteração é necessariamente diferente. Perde o `ctx`: usa `lssit%gpu`. |
| `cgm_kernels.cpp`, `hipsparse_generic.f90` | idem | Mantidos. O kernel de scatter pode sair, a depender de D4. |
| `cutestma.f90` | dois drivers | **Um só fonte**, linkado contra a biblioteca CPU e contra a GPU (`cutestma_cpu`/`cutestma_gpu`, D2(b)). |
| `gencan_gpu.f90`, `cutestma_gpu.f90`, `hess_csr.f90`, `cutest_iface.f90`, `gencan_gpu_context.f90` | — | **Removidos.** |
| `toyprob.f90` | fora do build | Volta ao build, como exemplo e como teste rápido dos dois modos. |

Saldo: 5 arquivos removidos e 1 novo (`lssit_gpu.f90`, que absorve o `hess_csr.f90`). As diferenças em relação ao original ficam concentradas em `lssit.f90` (estado e `preparehv`), `cgm.f90` (um desvio) e no mínimo de `gencan.f90`/`cutestma.f90`.

Fluxo, lado a lado:

```
original                          fase 6 (usegpu = .true.)
--------                          ------------------------
evalh      -> lssit%hrow/hcol/hval   (igual)
comphpre   (CPU)                     (igual)
preparehv  renumera + compacta       (igual) + espelha + coo2csr + envia CSR (lssit_gpu)
cgm        -> computehv / applyhpre  cgm -> cgm_gpu: vetores no device, SpMV, kernels de reducao
```

## 4. Decisões

### 4.0 O que foi decidido (22/11/2026)

| | Decisão | Observação |
|---|---|---|
| D1 | **(a)** pré-processador, opção `GENCAN_GPU` (`AUTO`/`ON`/`OFF`), `Fortran_PREPROCESS ON` | — |
| D2 | **(b)** escolha só em tempo de compilação, sem mudar a assinatura de `genunc`/`gencan` | **Provisória.** A opção (a) (argumento `usegpu`) depende de o orientador aceitar a mudança de assinatura **[ABERTO: perguntar ao orientador]**. O código deve facilitar a troca futura (ver abaixo). |
| D3 | estado da GPU dentro do `lssit_type` | — |
| D4 | espelhar + `coo2csr` in-place + enviar o CSR | **Medir** contra o scatter da fase 5 e procurar melhorias antes de fechar (ver abaixo). |
| D5 | código em `gencan-cg-standalone-gpu/` | — |
| D6 | branch `fase6`, merge em `main` após a validação final | — |

**Como fica D2(b) na prática.** Um só código-fonte e uma só árvore de build, compilando os fontes do Gencan duas vezes:
- biblioteca CPU: `gencan.f90`, `lssit.f90`, `cgm.f90` sem a macro `GENCAN_GPU`;
- biblioteca GPU (só se o HIP for encontrado): os mesmos fontes com `-DGENCAN_GPU`, mais `lssit_gpu.f90`, `cgm_gpu.f90` e os kernels;
- cada biblioteca com seu diretório de módulos (`Fortran_MODULE_DIRECTORY`), porque as duas geram `lssit.mod`, `modcgm.mod` etc.;
- o mesmo `cutestma.f90` (e `toyprob.f90`) linkado contra cada uma, gerando `cutestma_cpu` e `cutestma_gpu`. O `run-cutest-gpu.sh` quase não muda.

No modo GPU, `usegpu` é uma constante do `lssit` (`.true.` quando a macro está definida). Se o orientador aprovar D2(a), basta transformar essa constante em argumento de `genunc`/`gencan` e deixar uma biblioteca só.

Vantagem extra: o executável CPU sai dos mesmos fontes sem a macro, e é ele que se compara bit a bit com `gencan-cg-standalone/`.

**Como fica D4 na prática.** O passo 3 implementa o caminho proposto (espelhar + `coo2csr` + enviar o CSR) e mede, por passo de Newton, o tempo de `preparehv` (CPU), o envio e o SpMV, no `test_hess_csr` com `n` grande e em problemas CUTEst maiores. Alternativas a comparar se o custo pesar:
- (i) caminho proposto, CSR completo enviado a cada passo de Newton;
- (ii) com a face e o `hne` iguais, guardar a permutação que o `coo2csr` produziu e enviar só os `hne` valores do triângulo, espalhados na GPU como na fase 5;
- (iii) outras que aparecerem na medição.

A escolha final, com os números, vai para o changelog e para o PLANO §5.

### 4.1 Opções discutidas (proposta de 21/11/2026)

**D1. Como compilar sem HIP.**
- (a) **Recomendado:** pré-processador. Uma opção do CMake `GENCAN_GPU` (`AUTO`/`ON`/`OFF`), em que `AUTO` liga a GPU só se encontrar HIP/hipfort. Os arquivos continuam `.f90`, com `Fortran_PREPROCESS ON` no CMake (não precisa renomear para `.F90`). Os trechos de GPU ficam em `#ifdef GENCAN_GPU`. Sem HIP, os arquivos de GPU nem entram no build.
- (b) Alternativa: módulos "stub" com a mesma interface, escolhidos na hora de linkar. Os arquivos originais ficam sem `#ifdef`, mas surgem arquivos stub e a GPU fica menos visível na leitura.
- Motivo da recomendação: com (a), o `diff` contra a referência mostra exatamente onde a GPU entra, o que é o que o changelog precisa.

**D2. Como escolher CPU ou GPU ao rodar.**
- (a) **Recomendado:** um `logical :: usegpu` novo nos argumentos de `genunc`/`gencan`, no mesmo estilo de `truehp` e `precond`, definido no `cutestma`/`toyprob`. Se o build não tem GPU, `usegpu = .true.` vira CPU com um aviso. Se `truehp = .false.`, também vira CPU com aviso, porque a secante está fora do escopo por enquanto.
- Vantagem de (a): o mesmo executável roda CPU e GPU, então o `run-cutest-gpu.sh` compila uma vez por problema e roda duas. O `cutestma` pode aceitar `cpu`/`gpu` na linha de comando para os scripts.
- (b) Alternativa: só em tempo de compilação, sem mudar a interface do Gencan. Mas aí são necessários dois builds para comparar.

**D3. Onde fica o estado da GPU.**
- **Recomendado:** como componente do `lssit_type`, sob `#ifdef`. Isso atende ao "o `lssit` guarda a Hessiana" e elimina o `ctx` das assinaturas.
- Os vetores de trabalho do CG no device são alocados uma vez no `lssitini`, assim como hoje no `cgm_gpu_init`.

**D4. Adaptação do `coo2csr`.**
- Depois do `preparehv` original, o COO da face está em `lssit` (um triângulo, numeração `1..nfr`). O hipSPARSE precisa da matriz cheia, porque a API genérica do cuSPARSE não aceita "simétrica, só um triângulo" no SpMV, o que quebraria o NVIDIA.
- Proposta:
  - espelhar as entradas fora da diagonal no fim dos próprios vetores do `lssit`, que passam a ter capacidade `2*hnnzmax` quando `usegpu`;
  - chamar o `coo2csr` in-place, com `asta`/`alen` de tamanho `nfr`, que viram o `row_ptr`;
  - enviar o CSR.
- O padrão continua sendo reaproveitado: se a face e o `hne` não mudaram, a estrutura que sai do `coo2csr` é a mesma, porque a permutação depende só da sequência de linhas, que é fixa no CUTEst. Nesse caso só os valores sobem, e o descritor, o buffer e o preprocess continuam valendo, como hoje.
- Diferença de custo em relação à fase 5: hoje sobem `hne` doubles por passo de Newton e o scatter roda na GPU. Na proposta sobem cerca de `2*hne` doubles já em CSR, e o `coo2csr` roda na CPU (O(nnz), da mesma ordem do `preparehv`). Isso precisa ser medido (PLANO §7.2). Se pesar, a otimização de desempenho fica justificada e entra no changelog.
- A verificar: a ordem das colunas dentro de cada linha não é garantida. O rocSPARSE já funciona assim hoje (o `hess_csr` também não ordena). No cuSPARSE ainda não foi testado.

**D5. Onde fica o código.**
- **Recomendado:** continua em `gencan-cg-standalone-gpu/`, partindo da fase 5.
- `gencan-cg-standalone/` continua sendo a referência intocada e a base do changelog.
- O nome do diretório pode mudar no fim. Isso fica para depois.

**D6. Branch.**
- **Recomendado:** um branch `fase6`, com commits pequenos e merge em `main` quando a validação final passar. O `main` continua com a fase 5 validada como referência.

## 5. Passos (cada um termina com testes e um commit que atualiza o changelog)

1. **Preparação.** Criar o branch. Criar o `docs/CHANGELOG_GPU.md` com o formato combinado. Registrar a linha de base: as tabelas da fase 5 e os testes sintéticos.
2. **CMake.** Opção `GENCAN_GPU` (AUTO/ON/OFF), `Fortran_PREPROCESS`, as duas bibliotecas com diretórios de módulos separados (D2(b)) e os alvos `toyprob_{cpu,gpu}` e `cutestma_{cpu,gpu}` a partir de um só fonte cada. Conferir que `-DGENCAN_GPU=OFF` compila sem tocar em HIP.
   - Validação: o build OFF gera `toyprob`/`cutestma` só CPU, com resultados **bit a bit iguais** aos da referência.
3. **`lssit` + `lssit_gpu`.**
   - Tipo GPU dentro do `lssit_type`, com `lssitini`/`lssitend`.
   - Ramo GPU do `preparehv` com o espelhamento + `coo2csr` + envio.
   - SpMV, migrando o que hoje está em `hess_csr.f90`.
   - Validação: o `test_hess_csr` passa a testar `preparehv` com `usegpu` contra `computehv` da CPU, incluindo o caso "mesmo padrão, valores novos".
4. **`cgm` + `cgm_gpu`.** O desvio no `cgm`, e o `cgm_gpu` lendo `lssit%gpu` em vez do `ctx`.
   - Validação: `test_cgm_gpu`, 9/9 cenários.
5. **`gencan` + `cutestma`.**
   - `gencan.f90` sem `ctx`, com a GPU decidida pela macro (D2(b)).
   - Remover `gencan_gpu.f90` e `cutestma_gpu.f90`.
   - Conferir que o `run-cutest-gpu.sh` continua funcionando com os dois executáveis vindos do mesmo fonte.
6. **Limpeza.**
   - Remover o código morto e os scripts obsoletos.
   - Atualizar o `README.md` do diretório, o `CLAUDE.md` e o PLANO.
7. **Validação final.**
   - Testes sintéticos.
   - `pequena.txt` nos dois modos de `precond`, comparado com `resultados/table_*_fase5_*.txt`.
   - Build com `GENCAN_GPU=OFF`.
   - Se possível, build e `pequena.txt` no NVIDIA (RTX 5070 Ti).
   - Guardar as tabelas em `resultados/table_{cpu,gpu}_fase6_*.txt`.

## 6. Critério de aceitação

- **CPU:** idêntica à referência, bit a bit. O caminho CPU não pode mudar de comportamento.
- **GPU:** nos 6 problemas bem comportados de `pequena.txt`, `iter`/`fcnt`/`gcnt`/`hcnt` iguais aos da fase 5 e `f` dentro de 1e-6.
  - A ordem das entradas no CSR muda (agora vem do `coo2csr`), o que muda a ordem das somas no SpMV. O BQPGAUSS pode variar, como já documentado (PLANO §6, §8.10).
  - Se um problema bem comportado divergir, é preciso investigar antes de seguir.
- **Testes sintéticos:** 5/5 e 9/9 OK.
- **Build OFF:** compila numa máquina sem ROCm/CUDA.

## 7. Riscos

- `#ifdef` em excesso pode prejudicar a analogia. Regra: blocos curtos, só nos pontos de §3. Lógica de GPU com mais de algumas linhas vai para `lssit_gpu.f90`/`cgm_gpu.f90`.
- Mudar a assinatura de `genunc`/`gencan` (D2a) afetaria qualquer usuário externo do Gencan. Por isso D2(b) por enquanto, até a resposta do orientador.
- Compilar os mesmos fontes duas vezes (D2b) exige diretórios de módulos separados no CMake. Se misturar, um alvo pode pegar o `.mod` do outro sem erro aparente.
- Capacidade `2*hnnzmax` no `lssit` quando `usegpu`: mais memória na CPU (pequeno perto do resto).
- Desempenho do envio de CSR completo contra o scatter (D4): medir antes de decidir otimizar.
