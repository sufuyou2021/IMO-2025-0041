import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Path to Excel file
file_path = "C:/Users/mahler/Desktop/Figure S6.xlsx"

# Configuration for 4 target sheets and their corresponding group colors
sheet_configs = {

    "FS6_OFT_6W3": {
        "N_6W3": "#F6F6F6",
        "S_6W3": "#A0A0A0",
        "EA_6W3": "#3D3D3D"
    },
    "FS6_OFT_9W3": {
        "N_9W3": "#F6F6F6",
        "S_9W3": "#A0A0A0",
        "EA_9W3": "#3D3D3D"
    }
}

# Target column names to be used for plotting (standardized for comparison)
target_columns = [
    "Crossing Number",
    "Center Duration (sec)",
    "Center Distance (cm)",
    "Total Distance (cm)"
]

# Load the Excel file
xls = pd.ExcelFile(file_path)

# Process each sheet one by one
for sheet_name, group_colors in sheet_configs.items():
    # Read the sheet into a DataFrame
    df = pd.read_excel(xls, sheet_name=sheet_name)

    # Clean column names: strip spaces, remove symbols, replace space with underscore
    df.columns = df.columns.str.strip().str.replace(r'[^\w\s\(\)]', '', regex=True).str.replace(' ', '_')

    # Loop through each column (skip the first one which is assumed to be group column)
    for column in df.columns[1:]:
        # Only plot if the column matches target list
        if column.replace("_", " ") not in target_columns:
            continue

        # Create a figure for each boxplot
        fig, ax = plt.subplots(figsize=(8, 8))

        # Draw boxplot with black borders and custom fill colors
        sns.boxplot(
            x="gp", y=column, data=df,
            palette=group_colors,
            ax=ax, width=0.3, linewidth=6,
            boxprops=dict(edgecolor='black', facecolor='none', linewidth=6),
            whiskerprops=dict(color='black', linewidth=6),
            capprops=dict(color='black', linewidth=6),
            medianprops=dict(color='black', linewidth=6)
        )

        # Apply fill color to each box patch
        for patch, group in zip(ax.patches, df['gp'].unique()):
            patch.set_facecolor(group_colors.get(group, 'white'))
            patch.set_edgecolor("black")

        # Set title and axis labels
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel(column.replace("_", " "), fontsize=36, labelpad=20)

        # Adjust tick size and color
        ax.tick_params(axis='both', labelsize=30, width=6, colors="black")

        # Style plot borders
        for spine in ax.spines.values():
            spine.set_linewidth(6)
            spine.set_color("black")

        # Save plot as PDF
        output_name = f"{column}_{sheet_name}_boxplot_final.pdf"
        plt.tight_layout()
        plt.savefig(output_name, format="pdf", bbox_inches="tight")
        plt.close()
