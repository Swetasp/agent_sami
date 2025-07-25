# main.py
from agents.planner import Planner
from agents.executor import Executor
from agents.evaluator import Evaluator  # you can keep yours / stub it
from agents.memory import GlobalMemory
from code_sandbox import CodeSandbox
from utils.data_loader import load_csv_as_adata

# For now, use a local Ollama LLM. Replace with Groq/OpenAI if needed.
from langchain_community.llms import Ollama

def main():
    task = input("Enter your analysis task: ")
    data_path = input("Enter your data file path: ")

    # Load & summarize data
    adata = load_csv_as_adata(data_path)
    data_representation = str(adata)
    print(data_representation)
    # print("Data shape:", adata.shape)

    # data_signature = f"adata.shape={adata.shape}, obs={len(adata.obs.columns)}, var={len(adata.var.columns)}"

    # === Wire agents ===
    llm = Ollama(model="llama3", base_url="http://localhost:11434", temperature=0.0)
    memory = GlobalMemory()
    planner = Planner(llm, data_representation)
    sandbox = CodeSandbox("outputs/analysis.ipynb")
    executor = Executor(llm, data_path, memory, sandbox)
    evaluator = Evaluator(llm)

    # === Plan ===
    steps = planner.plan(task)
    print("\nPlan:")
    for i, s in enumerate(steps, 1):
        print(f"{i}. {s['description']}")

    # === Execute ===
    for i, step in enumerate(steps, 1):
        print(f"\n--- Executing step {i}: {step['description']} ---")
        result = executor.execute_step(step["description"])
        print("Success:", result.success)
        print("Stdout:\n", result.stdout)
        print("Stderr:\n", result.stderr)
        print("Notebook:", result.notebook_path or "not saved")

if __name__ == "__main__":
    main()
