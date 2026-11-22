#!/usr/bin/env python3
"""
Gera uma lista de problemas CUTEst filtrados pela classificação oficial.

Os filtros aceitam as strings descritivas usadas pelo pycutest:

  objetivo:
    none | constant | linear | quadratic | sum of squares | other

  restrições:
    unconstrained | fixed | bound | adjacency | linear | quadratic | other

Uso:
    python list_cutest_problems.py
    python list_cutest_problems.py --constraints unconstrained --max-n 1000
    python list_cutest_problems.py --objective "sum of squares" --constraints unconstrained
    python list_cutest_problems.py --constraints bound --min-n 10000 -o grandes_bound.txt
    python list_cutest_problems.py --help

Observação sobre n: na classificação, um problema de dimensão variável (n
escolhido no SIF) satisfaz qualquer filtro de n. Com --min-n, cada candidato é
decodificado com o sifdecoder (como no build) e o n do tamanho padrão é lido do
OUTSDIF.d; aí --min-n e --max-n valem para esse n. Requer SIFDECODE e MASTSIF
(source env_cutest.sh).
"""
import argparse
import os
import subprocess
import sys
import tempfile

try:
    import pycutest
except ImportError:
    print("pycutest não encontrado. Instale com: pip install pycutest", file=sys.stderr)
    sys.exit(1)


OBJECTIVES = ["none", "constant", "linear", "quadratic", "sum of squares", "other"]
CONSTRAINTS = ["unconstrained", "fixed", "bound", "adjacency",
               "linear", "quadratic", "other"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Lista problemas CUTEst por classificação (padrão: unconstrained).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--objective", default=None, choices=OBJECTIVES,
                   help="Tipo de objetivo (padrão: sem filtro)")
    p.add_argument("--constraints", default="unconstrained", choices=CONSTRAINTS,
                   help="Tipo de restrição (padrão: unconstrained)")
    p.add_argument("--max-n", type=int, default=None,
                   help="Número máximo de variáveis (n)")
    p.add_argument("--min-n", type=int, default=None,
                   help="Número mínimo de variáveis no tamanho padrão do SIF "
                        "(decodifica cada candidato; ver observação abaixo)")
    p.add_argument("--max-m", type=int, default=None,
                   help="Número máximo de restrições gerais (m)")
    p.add_argument("-o", "--output", default="problems.txt",
                   help="Arquivo de saída (padrão: problems.txt)")
    return p.parse_args()


def n_padrao(name):
    """n do tamanho padrao do SIF: decodifica num diretorio temporario e le o OUTSDIF.d."""
    sifdecoder = os.path.join(os.environ["SIFDECODE"], "bin", "sifdecoder")
    sif = os.path.join(os.environ["MASTSIF"], name + ".SIF")
    with tempfile.TemporaryDirectory() as d:
        try:
            r = subprocess.run([sifdecoder, sif], cwd=d, capture_output=True,
                               text=True, timeout=600)
        except subprocess.TimeoutExpired:
            return None
        out = os.path.join(d, "OUTSDIF.d")
        if r.returncode != 0 or not os.path.exists(out):
            return None
        with open(out) as f:
            for linha in f:              # primeiro inteiro do arquivo = n
                campos = linha.split()
                if campos:
                    try:
                        return int(campos[0])
                    except ValueError:
                        return None
        return None


def main():
    args = parse_args()

    kwargs = {"constraints": args.constraints}
    if args.objective is not None:
        kwargs["objective"] = args.objective
    if args.max_n is not None or args.min_n is not None:
        kwargs["n"] = (args.min_n or 0, args.max_n or 10**9)
    if args.max_m is not None:
        kwargs["m"] = (0, args.max_m)

    print(f"Filtrando CUTEst: {kwargs}...")
    names = sorted(pycutest.find_problems(**kwargs))

    if args.min_n is not None:
        for v in ("SIFDECODE", "MASTSIF"):
            if v not in os.environ:
                print(f"{v} nao definida: rode source env_cutest.sh", file=sys.stderr)
                sys.exit(1)
        print(f"Decodificando {len(names)} candidatos para ler o n padrao...")
        kept = []
        for name in names:
            n = n_padrao(name)
            ok = n is not None and n >= args.min_n and (args.max_n is None or n <= args.max_n)
            print(f"  {name:12s} n = {n if n is not None else 'falhou'}{'  <-' if ok else ''}")
            if ok:
                kept.append(name)
        names = kept

    with open(args.output, "w") as f:
        f.write("\n".join(names) + "\n")

    print(f"{len(names)} problemas encontrados; lista em {args.output}")


if __name__ == "__main__":
    main()