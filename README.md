# **SunLab Agent – Automated Multi-Omics Data Analysis Backend** 

SunLab Agent is an **LLM-driven backend framework** for automated biochemical and omics data analysis. It is inspired by **CellAgent** and designed to handle **multi-step analysis pipelines (normalization, clustering, integration, visualization)** through a **multi-agent architecture**.


## **Key Features**

- **Multi-Agent Architecture**:
  - **Planner** – Breaks down tasks into logical analysis steps.
  - **Executor** – Runs generated code and tools (normalization, clustering, etc.).
  - **Evaluator** – Evaluates the output and improves code when necessary.
  - **Global Memory** – Tracks and stores past steps and code executions.

- **Built-in Tools**:
  - `normalisation_tool` – Perform standard data normalization.
  - `clustering_tool` – Supports UMAP and Leiden clustering.
  - `integration_tool` – Align clusters from multiple datasets.

- **Data Compatibility**:
  - Supports **CSV**, **h5ad**, and other single-cell or multi-omics data formats.
  - Built on **Scanpy**, **Pandas**, and **Plotly** for powerful data analysis and visualization.

- **LLM Support**:
  - Works with **Groq, OpenAI, or local LLMs (e.g., Ollama)**.
  - Can dynamically plan, generate, and debug Python/R analysis code.

