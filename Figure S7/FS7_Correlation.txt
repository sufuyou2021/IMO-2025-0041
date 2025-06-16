import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr
import os

# Custom colormap: negative correlations = skyblue, positive = violet
custom_cmap = LinearSegmentedColormap.from_list(
    'skyblue_violet', ['skyblue', 'white', 'violet']
)

# Significance threshold for drawing bold borders
significance_threshold = 0.05

def generate_triangle_correlation_matrix_pdf(data, sheet_name, triangle='upper'):
    # Calculate Pearson correlation matrix
    corr_matrix = data.corr(method='pearson')

    # Initialize p-value matrix with 1s
    p_values_matrix = pd.DataFrame(np.ones_like(corr_matrix), index=corr_matrix.index, columns=corr_matrix.columns)

    # Calculate p-values for each correlation pair
    for i in range(len(corr_matrix.columns)):
        for j in range(i + 1, len(corr_matrix.columns)):
            col1 = corr_matrix.columns[i]
            col2 = corr_matrix.columns[j]
            _, p_value = pearsonr(data[col1].dropna(), data[col2].dropna())
            p_values_matrix.iloc[i, j] = p_value
            p_values_matrix.iloc[j, i] = p_value

    # Create mask: upper or lower triangle
    mask = np.tril(np.ones_like(corr_matrix, dtype=bool)) if triangle == 'upper' else np.triu(np.ones_like(corr_matrix, dtype=bool))

    # Draw heatmap
    plt.figure(figsize=(14, 12))
    ax = sns.heatmap(
        corr_matrix,
        mask=mask,
        annot=True,
        fmt='.2f',
        cmap=custom_cmap,
        cbar=True,
        vmin=-1, vmax=1,
        annot_kws={"size": 10, "color": "black"},
        linewidths=0,
        linecolor='gray'
    )

    # Format axis labels
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, fontstyle='italic', rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontstyle='italic', rotation=0)

    # Draw bold black borders around significant correlations
    for i in range(len(corr_matrix)):
        for j in range(len(corr_matrix)):
            if (triangle == 'upper' and j > i) or (triangle == 'lower' and j < i):
                if p_values_matrix.iloc[i, j] < significance_threshold:
                    ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=3))

    # Save plot as PDF
    plt.tight_layout()
    pdf_path = f'correlation_matrix_{sheet_name}_{triangle}_full_colorbar.pdf'
    plt.savefig(pdf_path, format='pdf')
    plt.close()
    print(f"✅ Saved: {pdf_path}")

# === Execution block ===
file_path = os.path.expanduser("~/Desktop/Figure S7.xlsx")

# Target sheets to process
sheets = [
    'FS7A_Correlation_APP6',
    'FS7B_Correlation_WT6',
    'FS7C_Correlation_APP9',
    'FS7D_Correlation_WT9'
]

# Generate both upper and lower triangle plots for each sheet
for sheet in sheets:
    data = pd.read_excel(file_path, sheet_name=sheet)
    generate_triangle_correlation_matrix_pdf(data, sheet, 'upper')
    generate_triangle_correlation_matrix_pdf(data, sheet, 'lower')
