import pandas as pd


def filter_entries(df):
    return df[
        (df["is_organic"] == True) &
        (df["has_metal"] == False) &
        (df["num_components"] == 1)
    ].copy()


def count_forms(df_filtered):
    df = df_filtered.dropna(subset=["inchi"]).copy()
    df["num_forms"] = df.groupby("inchi")["identifier"].transform("count")
    return df


def run_processing(input_path="csd_all.csv", output_path="csd_counted.csv"):
    df = pd.read_csv(input_path)

    df_filtered = filter_entries(df)
    print(f"Filtered: {len(df_filtered)} entries")

    df_counted = count_forms(df_filtered)
    print(f"Unique molecules: {df_counted['inchi'].nunique()}")

    df_counted.to_csv(output_path, index=False)
    return df_counted


if __name__ == "__main__":
    run_processing()
