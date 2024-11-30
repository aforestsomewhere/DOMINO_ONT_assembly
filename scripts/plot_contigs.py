import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load command-line arguments
taxonkit_file = sys.argv[1]
assembly_file = sys.argv[2]
output_file = sys.argv[3]
custom_colors = sys.argv[4].split(",")  # Split the comma-separated colors into a list

# Load taxonkit and assembly data
taxonkit_df = pd.read_csv(taxonkit_file, sep='\t', header=None, 
                          names=["level", "seq_name", "TaxID1", "TaxID2", "taxonomy"])
assembly_df = pd.read_csv(assembly_file, sep='\t')
assembly_df.rename(columns={"#seq_name": "seq_name"}, inplace=True)

# Process taxonomy column to keep only the last part after the last semicolon, handling NaN values
taxonkit_df["taxonomy"] = taxonkit_df["taxonomy"].apply(lambda x: x.split(";")[-1].strip() if isinstance(x, str) else "Unknown")

# Merge data on the "seq_name" column, which represents contig identifiers
merged_df = pd.merge(assembly_df, taxonkit_df, on="seq_name", how="inner")

# Filter the DataFrame to only include contigs with length > 10000
merged_df = merged_df[merged_df["length"] > 10000]

# Limit to the top 10 most frequent taxonomies and group others as "Other"
top_taxonomies = merged_df["taxonomy"].value_counts().index[:10]  # Get the top 10 unique taxonomies
merged_df["taxonomy"] = merged_df["taxonomy"].apply(lambda x: x if x in top_taxonomies else "Other")

# Sort by length in descending order for better visualization (optional)
merged_df = merged_df.sort_values(by="length", ascending=False)

# Generate a color map for the limited taxonomies (including "Other" if present)
unique_taxonomies = merged_df["taxonomy"].unique()[:10]  # Limit to at most 10 categories, including "Other"
color_map = dict(zip(unique_taxonomies, custom_colors[:len(unique_taxonomies)]))

# Assign colors based on taxonomy, replacing any NaNs with a default color
default_color = "gray"
colors = merged_df["taxonomy"].map(color_map).fillna(default_color)

# Verify there are no NaN values in the colors
if colors.isnull().any():
    colors = colors.fillna(default_color)

# Plot setup
plt.figure(figsize=(12, 8))

# Bar plot for contig lengths
bars = plt.bar(
    merged_df["seq_name"], 
    merged_df["length"], 
    color=colors  # Ensure colors are non-NaN
)

# Add borders for circular contigs
for bar, is_circular in zip(bars, merged_df["circ."]):
    if is_circular == "Y":
        bar.set_edgecolor("fuchsia")  # Add a fuchsia border for circular contigs
        bar.set_linewidth(3)          # Set the border width

# Custom legend for taxonomy colors and circular outline
legend_handles = [mpatches.Patch(color=color_map.get(tax, default_color), label=tax) for tax in unique_taxonomies]
circular_patch = mpatches.Patch(facecolor="none", edgecolor="fuchsia", linewidth=1, label="Circular Contig")
legend_handles.append(circular_patch)
plt.legend(handles=legend_handles, loc="best", title="Taxonomy and Circularity")

# Annotate each bar with the coverage information
for bar, cov in zip(bars, merged_df["cov."]):
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{cov}x', ha='center', va='bottom',fontsize=8)


# Labels and Title
plt.xlabel("Contigs")
plt.ylabel("Length (bp)")
plt.title("Contig Lengths, Coverage, and Circularity")

# Rotate x-axis labels and adjust layout
plt.xticks(rotation=45, fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout()

# Save the plot
plt.savefig(output_file, format='png', dpi=300)
