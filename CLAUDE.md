# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Projeto: Algencan/Gencan em GPU (HIP)

## Objetivo

Iniciação científica FAPESP (IME-USP, orientador Ernesto G. Birgin). Implementar e avaliar, em HIP, o gradiente conjugado com produtos Hessiana-vetor dentro do Algencan, em dois cenários: Hessiana explícita e quocientes incrementais (secante, `truehp = .false.`).

**Escopo atual:** só Hessiana explícita (`truehp = .true.`), no Gencan isolado (`gencan-cg-standalone-gpu/`). O projeto começou em 01/08/2026. Em 21/11/2026 estávamos no mês 4 de 12.

**Ordem de trabalho combinada (PLANO §7.1 e §7.5):**
1. **[FEITO 21/11/2026]** Terminar e validar a fase 5 na estrutura atual. Ela é a linha de base numérica (`gencan-cg-standalone-gpu/resultados/table_*_fase5_*.txt`).
2. Fase 6, a implementação provavelmente final, pedida pelo orientador:
   - código único CPU/GPU que compila mesmo sem HIP;
   - estrutura alinhada ao Gencan original, funcionando por analogia (mesmos módulos, rotinas e fluxo), com adaptações de desempenho justificadas;
   - um só `cutestma`, com o `lssit` guardando a Hessiana e o mínimo de arquivos novos;
   - COO → CSR adaptado de `docs/coo2csr_orientador.f90`;
   - cada mudança em relação ao Gencan de referência registrada em `docs/CHANGELOG_GPU.md`.

O plano da fase 6 está em `docs/FASE6.md`, aprovado em 22/11/2026: decisões no §4.0 e passos 1–7 no §5. O trabalho é feito no branch `fase6`. A escolha CPU/GPU é provisoriamente só de compilação (D2(b)), porque o argumento `usegpu` depende do orientador aprovar a mudança de assinatura. A rota do `coo2csr` (D4) precisa ser medida contra o scatter da fase 5.

Leia antes de começar:
- `docs/PLANO.md`: estado atual, decisões, testes, pendências, armadilhas e cronograma. É a fonte da verdade do andamento.
- `docs/Proposta.pdf`: a proposta FAPESP (objetivos, cronograma de 12 meses, critérios de validação).

## Stack e ambiente

- Fortran 90 (gfortran 16) + **hipfort**. Kernels em C++/HIP (`.cpp`) compilados com `hipcc`. Bindings `bind(C)` próprios para a API genérica do hipSPARSE (`hipsparse_generic.f90`).
- CMake ≥ 3.21 com presets `amd` (ROCm, `/opt/rocm`) e `nvidia` (backend CUDA). O mesmo código-fonte serve aos dois.
- Máquina de desenvolvimento: AMD RX 9070 XT (`gfx1201`). Há também uma máquina com RTX 5070 Ti, mas o backend NVIDIA ainda não foi testado com o código novo.
- Antes de qualquer build com CUTEst, carregue o ambiente. Sem isso o CMake aborta:
  `source <raiz>/env_cutest.sh` (define `CUTEST`, `SIFDECODE`, `MASTSIF`, `MYARCH=pc64.lnx.gfa`).
- Python (análises e gráficos): `.venv` + `requirements.txt`.

## Comandos (em `gencan-cg-standalone-gpu/`)

```bash
source ../env_cutest.sh
cmake --preset amd -DKERNEL_FLAGS="--offload-arch=gfx1201"
cmake --build --preset amd --target test_cgm_gpu && ./build/amd/test_cgm_gpu   # cgm_gpu x cgm (cenarios sinteticos)

# um problema CUTEst, CPU x GPU
cmake --preset amd -DPROBLEM=ROSENBR
cmake --build --preset amd --target cutestma_cpu cutestma_gpu
(cd build/amd && ./cutestma_cpu | tail -15 && ./cutestma_gpu | tail -15)

# varredura (retomavel): roda CPU e GPU por problema -> table_cpu.txt / table_gpu.txt
PRESET=amd ./run-cutest-gpu.sh pequena.txt
paste -d'\n' table_cpu.txt table_gpu.txt | less -S

# so CPU (sem HIP; sem CUTEst se o env nao estiver carregado)
cmake --preset cpu && cmake --build --preset cpu && ./build/cpu/toyprob_cpu
```

- GPU e CUTEst são opcionais no CMake:
  - `GENCAN_GPU=AUTO|ON|OFF`: o padrão sem preset é `AUTO`; os presets `amd`/`nvidia` usam `ON` e o preset `cpu` usa `OFF`.
  - Sem `env_cutest.sh`, os alvos `cutestma_*` não são criados.
  - Toda mudança no CMake precisa ser conferida nos dois backends (ROCm e CUDA). Nesta máquina não há CUDA: teste o nvidia pelo menos em AUTO/ON (fallback/erro) e registre o build CUDA real como **[NÃO VERIFICADO]** até rodar na máquina NVIDIA.

- **O problema CUTEst é fixado quando o CMake configura.** `-DPROBLEM=X` roda o `sifdecoder` e gera `ELFUN.f GROUP.f RANGE.f [EXTER.f]` em `build/<preset>/`. Para trocar de problema é preciso reconfigurar.
- `run-cutest-gpu.sh` pula os problemas listados em `done_gencan.txt` (apague a linha para repetir um problema). Os logs ficam em `build/<preset>/logs/<PROBLEMA>.log`, e é lá que aparece o motivo quando o script diz "falhou compilar". `pequena.txt` tem 7 problemas e `problems.txt` tem 475.
- Colunas da tabela: `Problem SC IERR CPU(s) n f gpsupn iter fcnt gcnt hcnt`. SC 0 significa convergiu, SC 99 significa que o tempo limite estourou.
- Opções do CMake para quando algo não é encontrado: `-DHIPFORT_ROOT=...` (ou `-DHIPFORT_MOD_DIR`/`-DHIPFORT_LIB`) e `-DCUDAToolkit_ROOT=...`. O hipfort precisa ter sido compilado com o mesmo gfortran.
- `./run-tests.sh` compila e roda `test_hess_csr` e `test_cgm_gpu`.
- Para repetir uma lista no `run-cutest-gpu.sh`, esvazie `done_gencan.txt` e deixe só o cabeçalho nas tabelas.
- `table_cpu.txt`/`table_gpu.txt` são a saída de trabalho do script. Resultados que valem como referência são copiados para `gencan-cg-standalone-gpu/resultados/` com o nome `table_{cpu,gpu}_<fase>[_<variante>].txt` e citados no PLANO §6.
- Divergência de `iter` entre CPU e GPU: antes de tratar como bug, rode o controle CPU contra CPU com a ordem das somas invertida (PLANO §6 e §8, caso BQPGAUSS).
- Scripts e alvos desatualizados:
  - `compile.sh` e `run-cutest.sh` deste diretório usam alvos que não existem com esses nomes (`toyprob`, `cutestma`; hoje há `toyprob_cpu` e `cutestma_{cpu,gpu}`).
  - O `README.md` cita arquivos fora do build (`scalar_kernels.cpp`, `gencan_gpu_context.f90`).

## Arquitetura (resumo; detalhes no PLANO §4)

O fluxo é o do Gencan original: `genunc`/`gencan` (gencan.f90) escolhe a face livre → `evalh` → `comphpre` → `preparehv` (lssit.f90) prepara a Hessiana da face → `cgm` (cgm.f90) resolve `H d = -g` nas variáveis livres, usando `computehv`/`applyhpre`.

Fase 6 (branch `fase6`): um só código. A GPU entra por analogia, em três pontos sob `#ifdef GENCAN_GPU`:
- **`lssit.f90`:** o `lssit_type` ganha `usegpu` (padrão `.true.` no build GPU) e `gpu` (`lssit_gpu.f90`). No modo GPU (`usegpu`, `truehp`, `compressed`), o `preparehv` chama `lssit_gpu_prepare` **no lugar** da compactação da CPU.
- **`lssit_gpu.f90`:** é a Hessiana da face na GPU.
  - Ele parte do COO do `evalh` e só remonta o padrão quando a estrutura do COO ou a face mudam: CSR com espelho, por contagem por linha, e `pos1/pos2`.
  - A cada passo, envia os `hne` valores e os espalha por kernel.
  - Guarda também os vetores de trabalho do CG.
  - Toda a memória da GPU é alocada uma vez.
- **`cgm.f90`:** se `lssit%gpu%active`, chama `cgm_gpu` (mesma assinatura) e retorna.
- **`cgm_gpu.f90`:** CG com a mesma ordem de testes do `cgm`, uma sincronização por iteração e trabalho especulativo. As reduções são kernels próprios (`cgm_kernels.cpp`), e o `comphpre` continua na CPU.
- **`gencan.f90`, `cutestma.f90` e `toyprob.f90` são iguais à referência.**
- **CMake:**
  - `gencan_cpu` = os fontes sem a macro;
  - `gencan_gpu` = os mesmos fontes com `-DGENCAN_GPU`, mais `lssit_gpu`/`cgm_gpu`;
  - `cutestma_{cpu,gpu}` e `toyprob_{cpu,gpu}` saem do mesmo fonte;
  - os testes comparam `usegpu = .false.` contra a GPU no mesmo executável.
- **Removidos em 22/11/2026:** `gencan_gpu.f90`, `cutestma_gpu.f90`, `cutest_iface.f90`, `gencan_gpu_context.f90` (estão no histórico). O `hess_csr.f90` só serve à comparação de tempos do `test_hess_csr`.

Armadilhas que já custaram tempo:
- Fortran não diferencia maiúsculas de minúsculas (`d_d` ≠ `d_dg`, nunca `d_D`).
- Os `add_executable(cutestma_*)` precisam vir depois do bloco do sifdecoder no CMake.
- Procedimento interno (depois do `contains` de um programa) não pode ter o próprio `contains`. Funções auxiliares de teste vão para o nível do programa.
- No build GPU, o `cgm` recebe o `lssit` como `intent(in)`. Qualquer estado de GPU que o CG precise já tem que estar alocado em `lssit%gpu` antes dele.
- A lista completa está no PLANO §3 e §8.

## Regras de trabalho

- Antes de editar, proponha o plano da tarefa e espere aprovação.
- Uma tarefa por vez, com commits pequenos e mensagens claras em português.
- Nunca altere a versão CPU de referência sem pedir: `gencan-cg-standalone/` inteiro.
- Desde a fase 6, `lssit.f90`, `cgm.f90` e `gencan.f90` em `gencan-cg-standalone-gpu/sources/` recebem a GPU **só sob `#ifdef GENCAN_GPU`**, em blocos curtos. Sem a macro, o comportamento tem que ser o da referência. Depois de cada mudança nesses arquivos, rode `./compare-cpu-ref.sh`, que compara a saída do build `cpu` (`toyprob_cpu` e `cutestma_cpu` na `pequena.txt`) com a referência compilada pelos scripts dela. Tem que dar 0 diferenças; só a linha de tempo é ignorada.
- Toda rotina ou kernel novo precisa ser testado contra o código CPU original antes de entrar no Gencan.
- Toda mudança numérica deve ser comparada com a CPU usando o mesmo `precond` dos dois lados. Critério da proposta (§6.3): `iter` idêntico, e `f`, gradiente e multiplicadores dentro de 1e-6. Critério do PLANO: `iter`, `fcnt`, `gcnt`, `hcnt` e `cginfo` idênticos. Se um problema bem condicionado divergir, trate como bug.
- Mudanças de desempenho: meça antes e depois com tempo de parede.
- Não rode experimentos longos sem avisar. Isso inclui varreduras do `problems.txt` e problemas grandes. `pequena.txt` e os testes sintéticos podem rodar sem aviso.
- Não altere o escopo, os objetivos ou o cronograma (PLANO e Proposta) sem consultar.
- Estilo: comentários no código e scripts em português sem acentos. Os documentos em `docs/` têm acentos.

## Registro

- Decisões técnicas vão em `docs/PLANO.md` §5 e resultados em `docs/PLANO.md` §6.
- A partir da fase 6, todo commit que mexe no código GPU atualiza também `docs/CHANGELOG_GPU.md` (arquivo/rotina, o que mudou em relação a `gencan-cg-standalone/`, por quê, como foi validado).
- `docs/Proposta.pdf` não vai para o git.
- Ao fim de cada sessão, atualize no PLANO a data, o estado dos testes (§6), as pendências (§7) e as armadilhas novas (§8).
- Marque o que foi escrito e ainda não foi rodado como **[NÃO VERIFICADO]**, e o que depende de decisão como **[ABERTO]**.

## Não fazer

- Não commitar binários, dados grandes ou credenciais. Diretórios `build*`, matrizes do SuiteSparse e `mastsif` ficam fora do controle de versão.
- Não remover código existente sem confirmação.
- Não editar `Initial-Benchmarks/` nem os CSVs de resultados dos benchmarks. São históricos.
