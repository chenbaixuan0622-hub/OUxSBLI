import os
import re
import glob

def parse_fortran_file(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # 1. 行継続(&)を結合して1行にする（パースを簡単にするため）
    merged_lines = []
    current_line = ""
    for line in lines:
        stripped = line.rstrip()
        # Fortranのコメント行は結合から除外
        if stripped.lstrip().startswith('!') and not stripped.lstrip().startswith('!>') and not stripped.lstrip().startswith('!<'):
            continue
            
        if stripped.endswith('&'):
            current_line += stripped[:-1] + " "
        else:
            current_line += stripped
            merged_lines.append(current_line)
            current_line = ""

    parsed_data = []
    current_block = None
    last_doc_comment = ""

    # 正規表現のコンパイル
    # サブルーチン/関数 (CUDAの attributes(...) 等にも対応)
    re_sub_func = re.compile(r'^\s*(?:attributes\s*\([^)]+\)\s+)?(?:pure\s+|elemental\s+|recursive\s+)?(subroutine|function|module)\s+(\w+)\s*(\(.*?\))?', re.IGNORECASE)
    # 変数宣言 (例: real(8), intent(in), device :: Q(5,nx) !< comment)
    re_var = re.compile(r'^\s*([^!]+?)\s*::\s*([^!]+?)(?:!<\s*(.+))?$', re.IGNORECASE)
    # ブロックの終了
    re_end = re.compile(r'^\s*end\s+(subroutine|function|module)', re.IGNORECASE)
    # 直前の説明コメント (!>)
    re_doc_before = re.compile(r'^\s*!>\s*(.+)')

    for line in merged_lines:
        # !> コメントの取得
        doc_match = re_doc_before.match(line)
        if doc_match:
            last_doc_comment = doc_match.group(1).strip()
            continue

        # サブルーチン/関数/モジュールの開始
        block_match = re_sub_func.match(line)
        if block_match:
            b_type = block_match.group(1).lower()
            b_name = block_match.group(2)
            b_args = block_match.group(3) if block_match.group(3) else ""
            
            current_block = {
                'type': b_type,
                'name': b_name,
                'args': b_args,
                'description': last_doc_comment,
                'variables': []
            }
            parsed_data.append(current_block)
            last_doc_comment = "" # リセット
            continue

        # ブロックの終了
        if re_end.match(line):
            current_block = None
            continue

        # 変数の取得（ブロック内部にいる時のみ）
        if current_block and current_block['type'] in ['subroutine', 'function']:
            var_match = re_var.match(line)
            if var_match and not line.lstrip().lower().startswith('use '):
                v_type_attrs = var_match.group(1).strip()
                v_names = var_match.group(2).strip()
                v_comment = var_match.group(3).strip() if var_match.group(3) else ""
                
                # 不要な空白を削除
                v_type_attrs = ' '.join(v_type_attrs.split())
                
                current_block['variables'].append({
                    'names': v_names,
                    'type_attrs': v_type_attrs,
                    'comment': v_comment
                })
            last_doc_comment = "" # 変数以外の行が来たらリセット

    return parsed_data

def generate_markdown_api(target_dir=".", output_file="api.md"):
    # 拡張子が .f90, .cuf, .f03 のものを検索
    extensions = ['*.f90', '*.cuf', '*.f03']
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(target_dir, '**', ext), recursive=True))

    if not files:
        print(f"指定されたディレクトリ ({target_dir}) にFortranファイルが見つかりません。")
        return

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("# Project API Index\n")
        f.write("> **AI Agent Instruction:** This is an auto-generated API reference for the CUDA Fortran project. Use this to understand the project structure, modules, subroutines, and their arguments without reading the full implementation details.\n\n")

        for filepath in files:
            data = parse_fortran_file(filepath)
            if not data:
                continue # 抽出可能なAPIがなければスキップ
            
            # 相対パスでファイル名を出力
            rel_path = os.path.relpath(filepath, target_dir)
            f.write(f"## File: `{rel_path}`\n\n")

            for block in data:
                f.write(f"### `{block['type']} {block['name']}{block['args']}`\n")
                if block['description']:
                    f.write(f"- **Description:** {block['description']}\n")
                
                if block['variables']:
                    f.write("- **Arguments & Variables:**\n")
                    for var in block['variables']:
                        comment_str = f" : {var['comment']}" if var['comment'] else ""
                        f.write(f"  - `{var['names']}` (`{var['type_attrs']}`){comment_str}\n")
                f.write("\n")
                
    print(f"AI向けのAPI仕様書を生成しました: {output_file}")

if __name__ == "__main__":
    # カレントディレクトリ以下の全ファイルを走査して api.md を作成
    generate_markdown_api(target_dir="../3D_solver/src", output_file="./docs/api.md")

