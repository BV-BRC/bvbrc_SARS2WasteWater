import pandas as pd
import os

# check the <sample_id>_freyja_variants.tsv exisits 
# edit the sample_time_metadata csv accordingly

df = pd.read_csv("sample_time_metadata.csv")
mask = []
for index, row in df.iterrows():
    sample_id = row["Sample"].split("_")[0]
    # print(os.path.join(path, "/home", "file.txt"))

    variant_file_path = os.path.join("output", sample_id, row["Sample"])
    print(variant_file_path)
    if os.path.exists(variant_file_path):
        print('file exisits')
        mask.append(True)
    else:
        mask.append(False)
filtered_df = df[mask]

# drop index 
filtered_df.to_csv("edited_sample_time_metadata.csv", index=False)