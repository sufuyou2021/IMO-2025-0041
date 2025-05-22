import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Excel 路徑
file_path = "C:/Users/mahler/Desktop/Figure S5.xlsx"

# 四個目標 sheet 與對應顏色
sheet_configs = {
    "FS5A_OFT_6": {
        "WT_6": "#F6F6F6",
        "APP/PS1_6": "#A0A0A0"
    },
    "FS5B_OFT_9": {
        "WT_9": "#F6F6F6",
        "APP/PS1_9": "#A0A0A0"
    },
    "FS5C_OFT_6W3": {
        "N_6W3": "#F6F6F6",
        "S_6W3": "#A0A0A0",
        "EA_6W3": "#3D3D3D"
    },
    "FS5D_OFT_9W3": {
        "N_9W3": "#F6F6F6",
        "S_9W3": "#A0A0A0",
        "EA_9W3": "#3D3D3D"
    }
}

# 欄位標題（統一比對格式）
target_columns = [
    "Crossing Number",
    "Center Duration (sec)",
    "Center Distance (cm)",
    "Total Distance (cm)"
]

# 開始處理每一個 Sheet
xls = pd.ExcelFile(file_path)

for sheet_name, group_colors in sheet_configs.items():
    df = pd.read_excel(xls, sheet_name=sheet_name)
    df.columns = df.columns.str.strip().str.replace(r'[^\w\s\(\)]', '', regex=True).str.replace(' ', '_')

    for column in df.columns[1:]:
        if column.replace("_", " ") not in target_columns:
            continue

        fig, ax = plt.subplots(figsize=(8, 8))

        sns.boxplot(
            x="gp", y=column, data=df,
            palette=group_colors,
            ax=ax, width=0.3, linewidth=6,
            boxprops=dict(edgecolor='black', facecolor='none', linewidth=6),
            whiskerprops=dict(color='black', linewidth=6),
            capprops=dict(color='black', linewidth=6),
            medianprops=dict(color='black', linewidth=6)
        )

        for patch, group in zip(ax.patches, df['gp'].unique()):
            patch.set_facecolor(group_colors.get(group, 'white'))
            patch.set_edgecolor("black")

        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel(column.replace("_", " "), fontsize=36, labelpad=20)
        ax.tick_params(axis='both', labelsize=30, width=6, colors="black")
        for spine in ax.spines.values():
            spine.set_linewidth(6)
            spine.set_color("black")

        output_name = f"{column}_{sheet_name}_boxplot_final.pdf"
        plt.tight_layout()
        plt.savefig(output_name, format="pdf", bbox_inches="tight")
        plt.close()
