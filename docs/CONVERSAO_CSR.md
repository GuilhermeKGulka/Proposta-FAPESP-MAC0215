# Hessiana da face na GPU: conversão COO → CSR (fase 6)

> Para explicar ao orientador o que o `preparehv` faz no modo GPU, por que está correto e por que é eficiente.
> Código: `gencan-cg-standalone-gpu/sources/lssit_gpu.f90` (`lssit_gpu_prepare`, `build_pattern`) e o kernel `k_scatter` em `cgm_kernels.cpp`.
> Escrito em 22/11/2026. Números medidos numa RX 9070 XT. Histórico das versões em `docs/CHANGELOG_GPU.md`.

## 1. O problema

A cada passo de Newton, o Gencan precisa de produtos `H_ff v`, em que `H_ff` é a submatriz da Hessiana nas variáveis livres (a face, `frvar`, com `nfr` variáveis).

- **O que o Gencan entrega:**
  - O `evalh` escreve no `lssit` a Hessiana em COO (`hrow`, `hcol`, `hval`, `hne` entradas), com índices `1..n` e **um triângulo**. Cada entrada fora da diagonal representa `H_ij` e `H_ji`. Entradas repetidas são somadas.
  - Esse é o contrato que o `computetruehv` da CPU já segue:
    ```
    hv(i) = hv(i) + hval(k) * v(j)
    if (j /= i) hv(j) = hv(j) + hval(k) * v(i)
    ```
- **O que a GPU precisa:** a matriz `nfr × nfr` em CSR (`row_ptr`, `col`, `val`), **cheia**. A API genérica do cuSPARSE não faz SpMV de matriz "simétrica, só um triângulo", e o mesmo código precisa servir ao ROCm e ao CUDA.

No Gencan original, quem prepara a Hessiana da face é o `preparehv`. Ele renumera as livres para `1..nfr` e compacta o COO, jogando fora as entradas que tocam variáveis presas. No modo GPU, o `preparehv` continua sendo esse ponto, mas faz o trabalho abaixo **no lugar** da compactação da CPU.

## 2. O algoritmo

O trabalho tem duas partes: uma cara, feita só quando o padrão muda, e uma barata, feita a cada passo.

### 2.1 Montagem do padrão (`build_pattern`): só quando a estrutura do COO ou a face mudam

1. **Renumerar** (a mesma regra do `preparehv`): `renum(i) = k` se `i` é a k-ésima variável livre, e `0` se está presa. Custo O(n).
2. **Contar as entradas por linha**, já com o espelho. Para cada entrada `e` do COO com `r = renum(hrow(e))` e `c = renum(hcol(e))`, as duas diferentes de 0, soma 1 em `alen(r)` e, se `r ≠ c`, também em `alen(c)`. Custo O(hne).
3. **Acumular o início de cada linha:** `row_ptr(1) = 1` e `row_ptr(i+1) = row_ptr(i) + alen(i)`. Custo O(nfr).
4. **Colocar as entradas fora do lugar, com um cursor por linha:** `cur(i) = row_ptr(i)`. Para cada entrada `e` mantida:
   - `j = cur(r)`, `col(j) = c`, `pos1(e) = j`, `cur(r) = j + 1`;
   - se `r ≠ c`: `j = cur(c)`, `col(j) = r`, `pos2(e) = j`, `cur(c) = j + 1`.

   As entradas descartadas ficam com `pos1 = pos2 = 0`. Custo O(hne).
5. **Enviar** `row_ptr`, `col`, `pos1` e `pos2` à GPU e recriar o descritor, o buffer e o preprocess do SpMV.

Os passos 2 e 3 são exatamente os dois primeiros laços do `coo2csr` do orientador. A diferença está no passo 4: o `coo2csr` faz uma permutação **in-place** seguindo ciclos, e aqui a colocação é **fora do lugar**, porque o resultado precisa ir para outro vetor de qualquer forma (ver §4).

**Exemplo** (n = 4, variável 2 presa, então `nfr = 3` e `renum = [1, 0, 2, 3]`):

```
COO (triangulo)        e:   1      2      3      4      5
                   (i,j):  (1,1)  (1,2)  (1,3)  (3,4)  (4,4)
                 (r,c) =   (1,1)  sai    (1,2)  (2,3)  (3,3)

passo 2: alen = [2, 2, 2]     (linha 1: (1,1),(1,2) | linha 2: (2,1),(2,3) | linha 3: (3,2),(3,3))
passo 3: row_ptr = [1, 3, 5, 7]
passo 4: e=1 -> col(1)=1, pos1=1
         e=3 -> col(2)=2, pos1=2 ; espelho col(3)=1, pos2=3
         e=4 -> col(4)=3, pos1=4 ; espelho col(5)=2, pos2=5
         e=5 -> col(6)=3, pos1=6
CSR: row_ptr = [1,3,5,7], col = [1,2 | 1,3 | 2,3]
```

### 2.2 A cada passo de Newton (`lssit_gpu_prepare`)

1. **Verificar se o padrão mudou:** comparar `hne`, `frvar`, `hrow` e `hcol` com os da última montagem. Se mudou, montar o padrão (§2.1).
2. **Enviar os valores:** sobem os `hne` valores do COO original (`8·hne` bytes).
3. **Espalhar na GPU** (`k_scatter`), com uma thread por entrada `e`: `val(pos1(e)) = hval(e)` e, se `pos2(e) > 0`, `val(pos2(e)) = hval(e)`.

O SpMV do hipSPARSE/cuSPARSE é então chamado pelo CG (`cgm_gpu`) sobre esse CSR.

## 3. Por que está correto

1. **Mesma matriz que a CPU usa.**
   - Cada entrada mantida `(r, c)` com valor `h` vira `h` na posição `(r, c)` do CSR e, se `r ≠ c`, também em `(c, r)`. É exatamente o que o `computetruehv` faz no produto.
   - Entradas repetidas no COO viram entradas repetidas no CSR, e o SpMV as soma, como a CPU.
   - Entradas que tocam variáveis presas são descartadas pela mesma regra do `preparehv` (`renum = 0`).
   - Logo, `CSR · v = H_ff v`, a matriz que o `computehv` usaria depois do `preparehv` original.
2. **Mesma numeração dos vetores.** A renumeração é a do `preparehv`, ou seja, a mesma ordem de `ifr = pack(1..n, frvar)` com que o `cgm` monta `g`, `x` e os limites compactados.
3. **O CSR é bem formado.**
   - A linha `i` recebe exatamente `alen(i)` entradas no passo 4: as mesmas condições do passo 2, contadas na mesma ordem.
   - O cursor de cada linha anda de `row_ptr(i)` até `row_ptr(i+1) − 1` sem sair dela.
   - Cada posição `1..nnz` recebe **exatamente uma** entrada (direta ou espelho). Assim, `pos1/pos2` formam uma correspondência um a um entre as entradas mantidas (e seus espelhos) e as posições do CSR.
4. **O scatter não tem corrida.** Como cada posição do CSR tem um único dono, as threads do `k_scatter` escrevem em posições distintas, e todas as posições são escritas a cada passo.
5. **Reaproveitar o padrão é válido.** `row_ptr`, `col` e `pos1/pos2` dependem só de `hne`, `hrow`, `hcol` e `frvar`. Se os quatro são iguais aos da última montagem, o padrão é o mesmo e só os valores precisam subir. A comparação é feita sem supor nada sobre o `evalh`, então vale também para quem não use o CUTEst.
6. **A ordem das colunas dentro da linha não importa** para o SpMV, e por isso não há ordenação. O rocSPARSE aceita colunas fora de ordem: é o que roda desde a fase 4. No cuSPARSE ainda falta conferir na máquina NVIDIA **[NÃO VERIFICADO]**.
7. **Evidência numérica:**
   - `test_hess_csr`: `H_ff v` da GPU contra `preparehv` + `computehv` da CPU, com erro relativo ≤ 2e-16. Os casos cobrem face cheia, face parcial, face vazia, "mesmo padrão com valores novos" e "mesma face com outra estrutura de COO".
   - `test_cgm_gpu`: o CG inteiro, 9/9 cenários com o mesmo `iter` e `cginfo` da CPU.
   - CUTEst `pequena.txt`: resultado idêntico ao da fase 5 em todos os campos.

## 4. Por que é eficiente

**Custo.** A montagem do padrão é O(n + hne): três passadas sequenciais pelo COO e uma pelo vetor de linhas. Cada passo comum é O(hne): enviar `8·hne` bytes e um kernel com uma thread por entrada.

**Memória.** Tudo é alocado **uma vez**, na primeira chamada, com a capacidade máxima `(n, hnnzmax)`. Quando o padrão muda, só o descritor e o buffer do SpMV são refeitos (exigência do hipSPARSE).

**Números** (teste sintético, `n = 300000`, 900000 entradas fora da diagonal, face com 80% das variáveis; `nnz` do CSR = 1,39 M):

| Por passo de Newton | Mesmo padrão | Face nova |
|---|---|---|
| Versão atual (contagem fora do lugar) | **1,65 ms** | **16 ms** |
| Fase 5 (`hess_csr`, mesma ideia, fora do `lssit`) | 1,05 ms | 14 ms |
| Versão com o `coo2csr` in-place (branch `fase6-coo2csr`) | 1,7 ms | 46 ms (32 ms só no `coo2csr`) |
| Primeira tentativa: compactar na CPU (`preparehv` original) + espelho + `coo2csr` + enviar o CSR inteiro a cada passo | 55 ms | — |
| `preparehv` original, só CPU (referência) | 23,6 ms | 24 ms |

Para comparar: o SpMV dessa matriz custa 0,065 ms.

**Por que o `coo2csr` in-place ficou mais lento aqui.**
- **Os dois primeiros laços são os mesmos.** A diferença está na colocação.
- **A permutação por ciclos tem uma cadeia de dependências.** Cada entrada levada ao seu lugar expulsa a que estava lá, e a próxima posição só é conhecida depois de ler essa entrada. Cada passo é um acesso aleatório que o processador não consegue antecipar, em três vetores (linha, coluna e a origem de cada entrada, necessária para obter `pos1/pos2`).
- **A colocação fora do lugar não tem esse problema.** Ela lê o COO em sequência, e as escritas em cada linha também são sequenciais.
- **A vantagem do in-place não se aplica aqui.** Ela é não precisar de memória extra, mas o resultado precisa de qualquer forma de `col`, `pos1` e `pos2` em vetores próprios (o COO do `lssit` não é alterado). Então o in-place não economiza nada.

**Por que não compactar na CPU a cada passo.** Na primeira versão, o `preparehv` original rodava antes da conversão. Só essa compactação custa 24 ms por passo, e a GPU não precisa dela: com `pos1/pos2`, o filtro da face já fica embutido no padrão.

**Quanto o padrão muda na prática.** No Gencan, com a `pequena.txt`:
- sem caixa, o padrão é montado **uma vez**;
- com caixa, a face muda em quase todo passo de Newton (BQPGAUSS 109 de 111, TORSION1 23 de 23, OBSTCLAE 21 de 21).

Por isso o custo da montagem (16 contra 46 ms) importa nos problemas com caixa. A medição em problemas grandes está no PLANO §6.

## 5. Ponto em aberto: verificar se a estrutura mudou

A comparação do passo 2.2.1 lê `hrow`, `hcol` e os vetores guardados: cerca de 19 MB no exemplo, ~0,6 ms por passo, limitado pela banda de memória.
- **Comparar a face é barato:** `n` lógicos.
- **O caro é a estrutura do COO.** Só quem escreveu o `evalh` sabe se ela é fixa (no CUTEst é). A interface atual do Gencan não tem como receber essa informação.
- **Opções** (detalhes no `docs/FASE6.md` §4.0):
  - (1) manter a comparação, que está correta para qualquer `evalh`;
  - (2) supor estrutura fixa, como a fase 5;
  - (3) um argumento novo em `genunc`/`gencan` (ou no `evalh`) para o usuário avisar que a estrutura é fixa. Isso muda a assinatura, junto com o `usegpu` do D2(a).

Hoje vale a opção (1).
