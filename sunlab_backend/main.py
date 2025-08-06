# main.py
from agents.planner import Planner
from agents.executor import Executor
from agents.evaluator import Evaluator  # you can keep yours / stub it
from agents.memory import GlobalMemory
from code_sandbox import CodeSandbox
from utils.data_loader import load_csv_as_adata
from dotenv import load_dotenv
import os
import pandas as pd
import json

# For now, use a local Ollama LLM. Replace with Groq/OpenAI if needed.
from langchain_community.llms import Ollama

# from langchain_community.chat_models import ChatOpenAI
from langchain_openai import ChatOpenAI

def main():
    task = input("Enter your analysis task: ")
    data_path = input("Enter your data file path: ")

    df = pd.read_csv(data_path)

    summary = {
        "num_rows": df.shape[0],
        "num_columns": df.shape[1],
        "columns": df.columns.tolist(),
        "dtypes": df.dtypes.astype(str).to_dict()
    }

    column_summaries = {}

    for col in df.columns:
        if df[col].dtype in ['float64', 'int64']:
            column_summaries[col] = df[col].describe().to_dict()
        else:
            column_summaries[col] = {
                "type": str(df[col].dtype),
                "unique_values": df[col].nunique(),
                "most_common": df[col].value_counts().head(3).to_dict()
            }
    preview = df.head(3).to_dict(orient="records")

    data_summary = {
    "shape": df.shape,
    "column_types": summary["dtypes"],
    "column_summaries": column_summaries,
    "preview": preview
    }

    summary_string = json.dumps(data_summary, indent=2)
    # print(summary_string)

    # Load & summarize data
    adata = load_csv_as_adata(data_path)
    data_representation = str(adata)
    # print(data_representation)
    # print("Data shape:", adata.shape)

    # data_signature = f"adata.shape={adata.shape}, obs={len(adata.obs.columns)}, var={len(adata.var.columns)}"

    load_dotenv()
    JUPYTER_NOTEBOOK_PATH = os.environ.get("JUPYTER_NOTEBOOK_PATH")
    OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
    
    llm = ChatOpenAI(
    model="gpt-4o",
    temperature=0.0,
    openai_api_key=OPENAI_API_KEY,
    openai_api_base="https://api.ai.it.ufl.edu/v1",
    )
    global_memory = GlobalMemory()
    planner = Planner(llm, data_representation= summary_string)
    sandbox = CodeSandbox("outputs/analysis.ipynb")
    executor = Executor(llm, data_path, global_memory, sandbox)
    evaluator = Evaluator(llm)


    steps = planner.plan(task)
    print("\nPlan:")
    for i, s in enumerate(steps, 1):
        print(f"{i}. {s['description']}")

    for i, step in enumerate(steps, 1):
        print(f"\n--- Executing step {i}: {step['description']} ---")
        # result = executor.execute_step(step["description"])
        result = executor.execute_step(
            step_description=step["description"],
            user_requirements=task,
            data_description=data_representation,
            global_memory=global_memory,
        )
    
        global_memory.add_code(result.code)
        print("Success:", result.success)
        print("Stdout:\n", result.stdout)
        print("Stderr:\n", result.stderr)
        print("Notebook:", result.notebook_path or "not saved")
    executor.finalize_notebook()

    # Run evaluation
    evaluation = evaluator.evaluate_step(
        step_description=step["description"],
        code=result.code,
        stdout=result.stdout,
        stderr=result.stderr,
    )

    print("\nEvaluation:")
    print(f"- Passed: {evaluation.passed}")
    print(f"- Score: {evaluation.score}")
    print(f"- Explanation: {evaluation.explanation}")

if __name__ == "__main__":
    main()