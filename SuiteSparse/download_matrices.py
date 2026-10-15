#!/usr/bin/env python3
"""Baixa da SuiteSparse (formato Matrix Market) as matrizes listadas em um arquivo.

Uso:
    python download_matrices.py [lista.txt] [--dest DIR] [--max-nnz N] [--dry-run]

A lista tem um nome por linha ('#' inicia comentario); aceita 'grupo/nome' para
desambiguar (o grupo so serve para achar a matriz certa; a pasta final usa so o nome).

Layout gerado (o mesmo que os run_benchmark*.sh esperam):
    <dest>/MM/<nome>/<arquivos...>

O ssgetpy baixa e extrai em <dest>/<nome>/; este script move para <dest>/MM/<nome>/.

Uma matriz conta como "ja baixada" se existir <dest>/MM/<nome>/ com pelo menos um
arquivo dentro. Pastas deixadas em <dest>/<nome>/ por execucoes antigas sao movidas
para <dest>/MM/ (e contam como "movidas da raiz" no resumo). Pode rodar de novo a
vontade.

Antes de comecar, o script checa se o servidor da SuiteSparse esta servindo os
arquivos: ele as vezes fica offline e devolve uma pagina HTML no lugar do .tar.gz
(por isso os erros "Not a gzipped file (b'<h')"). Nesse caso o script aborta com
instrucoes, em vez de tentar baixar e falhar em todas as 238. Use
--skip-offline-check para pular essa checagem.
"""
import argparse
import os
import shutil
import sys
import urllib.request

import ssgetpy

HERE = os.path.dirname(os.path.abspath(__file__))

# Frases que aparecem na pagina quando o servico esta com o download desabilitado.
OFFLINE_MARKERS = (
    "THE COLLECTION IS OFFLINE",
    "COLLECTION IS OFFLINE",
)


def read_list(path):
    items = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                items.append(line)
    return items


def check_offline(timeout=15):
    """Devolve (offline, detalhe).

    offline = True  -> o site esta servindo o alerta de offline
    offline = False -> o site respondeu e nao parece estar com o alerta
    offline = None  -> nao foi possivel checar (sem rede, DNS, TLS, proxy, ...)
    """
    url = "https://sparse.tamu.edu/"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            head = resp.read(16000).decode("utf-8", errors="replace").upper()
    except Exception as e:
        return None, f"nao foi possivel contatar {url}: {e}"
    for marker in OFFLINE_MARKERS:
        if marker in head:
            return True, f"{url} esta exibindo o alerta de offline"
    return False, f"{url} respondeu normalmente"


def has_files(path):
    """True se 'path' e um diretorio com pelo menos um arquivo dentro (recursivo)."""
    for _root, _dnames, fnames in os.walk(path):
        if fnames:
            return True
    return False


def move_into_mm(name, dest, dry_run, tag=""):
    """Move <dest>/<name>/ para <dest>/MM/<name>/ se necessario.

    Devolve True se moveu (ou moveria, em dry-run). Se nao havia nada para mover
    devolve False. Se MM/<name>/ ja existe, apenas apaga a sobra em <dest>/<name>/.
    """
    mm_root = os.path.join(dest, "MM")
    src = os.path.join(dest, name)
    dst = os.path.join(mm_root, name)
    if not (os.path.isdir(src) and has_files(src)):
        return False
    if os.path.isdir(dst):
        # ja esta no lugar certo; a copia na raiz e sobra, apaga.
        print(f"{tag} [limpando] sobra em {src} (MM/{name} ja existe)")
        if not dry_run:
            shutil.rmtree(src, ignore_errors=True)
        return False
    print(f"{tag} [movendo] {src} -> {dst}")
    if not dry_run:
        os.makedirs(mm_root, exist_ok=True)
        shutil.move(src, dst)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lista", nargs="?", default=os.path.join(HERE, "matrices.txt"))
    ap.add_argument("--dest", default=HERE, help="pasta base (padrao: a deste script)")
    ap.add_argument("--max-nnz", type=int, default=None, help="pula matrizes com mais nnz que isso")
    ap.add_argument("--dry-run", action="store_true", help="so mostra o que seria baixado/movido")
    ap.add_argument("--skip-offline-check", action="store_true",
                    help="nao checa se o servidor esta offline antes de comecar")
    args = ap.parse_args()

    names = read_list(args.lista)
    mm_root = os.path.join(args.dest, "MM")
    print(f"{len(names)} matrizes na lista '{args.lista}' -> {mm_root}")

    if not args.skip_offline_check:
        print("Checando o servidor da SuiteSparse...")
        offline, detail = check_offline()
        if offline is True:
            print(f"[offline] {detail}")
            print("  O servico de download da SuiteSparse esta desabilitado no momento.")
            print("  Veja https://sparse.tamu.edu/ (o alerta fica no topo da pagina).")
            print("  Rode de novo quando o servico voltar, ou use --skip-offline-check para tentar mesmo assim.")
            return 2
        elif offline is False:
            print(f"[ok] {detail}")
        else:
            print(f"[aviso] {detail}; tentando de qualquer forma")

    done, skipped, adopted, too_big, missing, failed = 0, 0, 0, [], [], []
    interrupted = False
    try:
        for idx, item in enumerate(names, 1):
            group, _, name = item.rpartition("/")
            tag = f"[{idx}/{len(names)}]"
            dst = os.path.join(mm_root, name)

            # 1) ja existe em MM/<nome>/ com arquivos dentro?
            if os.path.isdir(dst) and has_files(dst):
                skipped += 1
                continue

            # 2) sobrou em <dest>/<nome>/ de uma execucao antiga? move para MM/.
            if move_into_mm(name, args.dest, args.dry_run, tag=tag):
                adopted += 1
                continue

            # 3) precisa baixar.
            found = [m for m in ssgetpy.search(name=name, limit=1000)
                     if m.name == name and (not group or m.group == group)]
            if not found:
                missing.append(item)
                print(f"{tag} [nao encontrada] {item}")
                continue
            m = found[0]
            if len(found) > 1:
                print(f"[aviso] '{name}' existe em varios grupos "
                      f"({', '.join(x.group for x in found)}); usando '{m.group}'")
            nnz = getattr(m, "nnz", None)
            if args.max_nnz and nnz and nnz > args.max_nnz:
                too_big.append(f"{m.group}/{m.name} (nnz={nnz})")
                continue

            print(f"{tag} [{'dry-run' if args.dry_run else 'baixando'}] "
                  f"{m.group}/{m.name}  nnz={nnz}", flush=True)
            if args.dry_run:
                continue
            try:
                m.download(format="MM", destpath=args.dest, extract=True)
                # ssgetpy extrai em <dest>/<name>/; move para <dest>/MM/<name>/.
                if not move_into_mm(name, args.dest, dry_run=False, tag=tag):
                    if os.path.isdir(dst) and has_files(dst):
                        pass  # ja caiu no lugar certo (versao diferente do ssgetpy)
                    else:
                        raise RuntimeError(
                            f"download nao deixou arquivos em "
                            f"{os.path.join(args.dest, name)} nem em {dst}"
                        )
                done += 1
            except Exception as e:
                failed.append(item)
                print(f"[falhou] {item}: {e}")
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrompido. Rode de novo para continuar de onde parou.")

    print("-" * 50)
    print(f"baixadas: {done} | ja existiam: {skipped} | movidas da raiz: {adopted} | "
          f"puladas por --max-nnz: {len(too_big)} | nao encontradas: {len(missing)} | "
          f"falhas: {len(failed)}")
    for label, lst in (("puladas (grandes)", too_big),
                       ("nao encontradas", missing),
                       ("falhas", failed)):
        if lst:
            print(f"{label}: {', '.join(lst)}")
    if interrupted:
        return 130
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())