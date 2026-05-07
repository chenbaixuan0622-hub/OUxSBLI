import pandas as pd
import json
import re
import os

def parse_numeric_advanced(val_str):
    if pd.isna(val_str): return 0.0
    s = str(val_str).strip()
    match = re.search(r'\d+\.?\d*', s)
    if match:
        return float(match.group(0))
    return 0.0

# 抽出したいターゲットセクションの定義（Nsight Computeのバージョンによる表記揺れを吸収）
TARGET_SECTIONS_MAP = {
    "speed of light": "GPU Speed Of Light Throughput",
    "speedoflight": "GPU Speed Of Light Throughput",
    "memory workload": "Memory Workload Analysis",
    "memoryworkload": "Memory Workload Analysis",
    "scheduler": "Scheduler Statistics",
    "warp state": "Warp State Statistics",
    "warpstate": "Warp State Statistics",
    "launch": "Launch Statistics",
    "occupancy": "Occupancy"
}

def get_target_section(section_name):
    """CSVのセクション名が抽出対象のセクションに該当するか判定する"""
    if pd.isna(section_name): return None
    s_lower = str(section_name).lower()
    for key, display_name in TARGET_SECTIONS_MAP.items():
        if key in s_lower:
            return display_name
    return None

def generate_metrics_focused_report():
    # 1. ターミナルから入力ファイル名を取得
    input_file = input("Nsight Computeから出力したCSVファイル名を入力してください（例: report.csv）: ").strip()
    
    if not os.path.exists(input_file):
        print(f"エラー: ファイル '{input_file}' が見つかりません。")
        return

    # 2. 出力ファイル名の生成 (拡張子を .json に変更)
    base_name = os.path.splitext(input_file)[0]
    output_json = base_name + ".json"

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"CSVの読み込みエラー: {e}")
        return

    df.columns = [c.strip() for c in df.columns]

    # カラム名の検知
    metric_col = next((c for c in df.columns if 'metric name' in c.lower()), 'Metric Name')
    val_col = next((c for c in df.columns if 'metric value' in c.lower()), 'Metric Value')
    unit_col = next((c for c in df.columns if 'metric unit' in c.lower()), 'Metric Unit')
    section_col = next((c for c in df.columns if 'section name' in c.lower()), 'Section Name')

    kernels_data = {}

    # 3. データの収集
    for _, row in df.iterrows():
        k_id = row.get('ID')
        if pd.isna(k_id): continue
        
        if k_id not in kernels_data:
            kernels_data[k_id] = {
                "kernel_id": int(k_id) if isinstance(k_id, (int, float)) else k_id,
                "kernel_name": row.get('Kernel Name', 'Unknown'),
                "duration_ms": 0.0,
                "metrics": {} # 指定されたセクションのメトリクスをここに格納
            }
            
        m_name = str(row.get(metric_col, '')).strip()
        m_val = str(row.get(val_col, '')).strip()
        m_unit = str(row.get(unit_col, '')).strip()
        m_name_lower = m_name.lower()

        # --- 実行時間 (Duration) の取得 ---
        if m_name_lower == 'duration' or ('duration' in m_name_lower and 'time' in m_name_lower):
            val = parse_numeric_advanced(m_val)
            unit_lower = m_unit.lower()
            
            if unit_lower in ['ns', 'nsecond']: val /= 1_000_000
            elif unit_lower in ['us', 'usecond']: val /= 1_000
            elif unit_lower in ['s', 'second']: val *= 1000
            
            kernels_data[k_id]['duration_ms'] = max(kernels_data[k_id]['duration_ms'], val)

        # --- 指定セクションのメトリクス取得 ---
        section_name = row.get(section_col)
        target_section = get_target_section(section_name)
        
        if target_section:
            if target_section not in kernels_data[k_id]['metrics']:
                kernels_data[k_id]['metrics'][target_section] = {}
            
            # 値と単位を結合して見やすくする (例: "45.2 %", "1024 byte")
            display_val = m_val
            if m_unit and m_unit.lower() not in ['none', 'n/a', 'null', '']:
                display_val = f"{m_val} {m_unit}"
                
            kernels_data[k_id]['metrics'][target_section][m_name] = display_val

    # 4. 実行時間が長い順にソートして抽出
    kernel_list = [k for k in kernels_data.values() if k['duration_ms'] > 0]
    kernel_list.sort(key=lambda x: x['duration_ms'], reverse=True)
    
    # 抽出するカーネルの数（必要に応じて変更してください）
    top_k = 5
    selected_kernels = kernel_list[:top_k]

    ai_context = {
        "source_file": input_file,
        "bottleneck_kernels": []
    }

    for k in selected_kernels:
        ai_context["bottleneck_kernels"].append({
            "kernel_id": k["kernel_id"],
            "kernel_name": k["kernel_name"],
            "duration_ms": round(k["duration_ms"], 6),
            "metrics": k["metrics"]
        })

    # 5. JSONとして出力
    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(ai_context, f, indent=2, ensure_ascii=False)
        
    print("-" * 40)
    print(f"解析成功:")
    print(f"  入力: {input_file}")
    print(f"  出力: {output_json}")
    print(f"  抽出セクション: {list(TARGET_SECTIONS_MAP.values())[:6]} ...")
    print("-" * 40)
    print("【抽出された重いカーネル Top】")
    for i, k in enumerate(ai_context["bottleneck_kernels"]):
        print(f"{i+1:2d}. {k['duration_ms']:10.6f} ms : {k['kernel_name']}")

if __name__ == "__main__":
    generate_metrics_focused_report()
