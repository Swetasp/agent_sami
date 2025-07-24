import pandas as pd
import scanpy as sc

def load_csv_as_adata(
    path: str,
    index_col: int | str = 0,
    store_non_numeric_in_obs: bool = True,
) -> sc.AnnData:
    """
    Load a generic CSV and convert it to AnnData.

    - First column is assumed to be the row index.
    - If the index contains duplicates, we make it unique by appending _0, _1, ...
    - Numeric columns -> adata.X
    - Non-numeric columns (optional) -> adata.obs
    """
    df = pd.read_csv(path)

    # Set the index to the first column (or a named column)
    if isinstance(index_col, int):
        df = df.set_index(df.columns[index_col])
    else:
        df = df.set_index(index_col)

    # ---- make index unique if needed ----
    if not df.index.is_unique:
        dup_counts = df.groupby(level=0).cumcount()
        df.index = df.index.astype(str) + "_" + dup_counts.astype(str)

    # Split numeric / non-numeric
    numeric_df = df.select_dtypes(include=["number"])
    non_numeric_df = df.select_dtypes(exclude=["number"])

    if numeric_df.empty:
        raise ValueError(
            "No numeric columns found in the CSV. "
            "Please check the file or adjust how you read it."
        )

    # Build AnnData
    adata = sc.AnnData(numeric_df.to_numpy())
    adata.obs_names = numeric_df.index.astype(str)
    adata.var_names = numeric_df.columns.astype(str)

    # Attach non-numeric columns as obs metadata if requested
    if store_non_numeric_in_obs and not non_numeric_df.empty:
        # Align on the (now unique) index
        non_numeric_df = non_numeric_df.reindex(numeric_df.index)
        adata.obs = pd.concat([adata.obs, non_numeric_df], axis=1)

    return adata
