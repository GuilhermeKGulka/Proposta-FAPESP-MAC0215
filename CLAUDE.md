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

O fluxo é o do Gencan original: `genunc`/`gencan` (gencan.f90) escolhe a face livre → `cgm` (cgm.f90) resolve `H d = -g` nas variáveis livres → `computehv`/`applyhpre` (lssit.f90). A versão GPU troca só o CG:

- `gencan_gpu.f90` é uma cópia do `gencan.f90` com mudanças mínimas (PLANO §4.4). No lugar de `preparehv` + `cgm`, chama `hess_csr_update` + `cgm_gpu`. O contexto GPU é local a `genunc`/`gencan`.
- `hess_csr.f90` monta, na GPU, o CSR `nfr×nfr` a partir do COO do `evalh` e do `frvar`, com a mesma renumeração de `preparehv` (que não é chamado). O padrão, o descritor e o buffer só são refeitos quando a face muda. Fora isso, apenas os valores são espalhados por kernel.
- `cgm_gpu.f90` tem a mesma interface e a mesma ordem de testes que `cgm`. Os vetores ficam no device e as reduções são kernels próprios em `cgm_kernels.cpp`. Na fase 5 há uma sincronização por iteração, com o trabalho da iteração seguinte lançado de forma especulativa. `comphpre` continua na CPU.
- CMake: `gencan_cpu` (referência intocada) ← `gpucore` (hess_csr + cgm_gpu) ← `gencan_gpu`. `kernels` é o objeto gerado pelo hipcc. `cutestma.f90` e `cutestma_gpu.f90` são os drivers CUTEst e precisam permanecer equivalentes, inclusive no valor de `precond`.

Armadilhas que já custaram tempo:
- Fortran não diferencia maiúsculas de minúsculas (`d_d` ≠ `d_dg`, nunca `d_D`).
- Os `add_executable(cutestma_*)` precisam vir depois do bloco do sifdecoder no CMake.
- Testes sintéticos precisam de `hess_csr_invalidate` antes de cada `hess_csr_update`.
- A lista completa está no PLANO §3 e §8.

## Regras de trabalho

- Antes de editar, proponha o plano da tarefa e espere aprovação.
- Uma tarefa por vez, com commits pequenos e mensagens claras em português.
- Nunca altere a versão CPU de referência sem pedir: `gencan-cg-standalone/` inteiro.
- Desde a fase 6, `lssit.f90`, `cgm.f90` e `gencan.f90` em `gencan-cg-standalone-gpu/sources/` recebem a GPU **só sob `#ifdef GENCAN_GPU`**, em blocos curtos. Sem a macro, o comportamento tem que ser o da referência. Depois de cada mudança nesses arquivos, compare a saída do build `cpu` (`toyprob_cpu` e `cutestma_cpu` na `pequena.txt`) com a referência compilada pelos scripts dela: tem que ser idêntica, fora a linha de tempo.
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
