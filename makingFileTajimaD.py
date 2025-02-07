#this paython code is for combining Tajima'sD calucluation into single file to analyze in R. If your output file contains nana, please remove nana by grep -v nan file before running this code.
import pandas as pd
import os
import matplotlib.pyplot as plt

# Define file paths and corresponding group names
files = {
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_bif_10K.Tajima.D_ed': 'Abifurcata',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_cyt_10K.Tajima.D_ed': 'Acytherea',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_hya_10K.Tajima.D_ed': 'Ahyacinthus',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_sub_10K.Tajima.D_ed': 'Asubulata',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_gem_10K.Tajima.D_ed': 'Agemmifera',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_hum_10K.Tajima.D_ed': 'Ahumilis',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_mon_10K.Tajima.D_ed': 'Amonticulosa',
    '/path/GemHumMonHyaRefGem_geno005_maf001_hwe_hyb_10K.Tajima.D_ed': 'Hybrids'
}

# Collect file paths and group names
data_frames = []
for file, group_name in files.items():
    if os.path.isfile(file):
        try:
            df = pd.read_csv(file, sep='\t')
            df['Group'] = group_name  # Add group name column
            data_frames.append(df[['CHROM', 'BIN_START', 'N_SNPS', 'TajimaD', 'Group']])
        except Exception as e:
            print(f"Error reading {file}: {e}")
    else:
        print(f"File does not exist: {file}")

if data_frames:
    # Combine all data into a single DataFrame
    combined_df = pd.concat(data_frames, ignore_index=True)

    # Display the first few rows of the combined data
    print(combined_df.head())

    # Save the combined data to a file (optional)
    combined_df.to_csv('combined_tajimaD_with_groups_gem_hya_ref_Agem_windo10K.csv', index=False)
