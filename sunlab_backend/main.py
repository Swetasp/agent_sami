from agents.planner import Planner
from agents.executor import Executor
from agents.evaluator import Evaluator
from agents.memory import GlobalMemory
from utils.data_loader import load_csv_as_adata
from langchain_community.llms import Ollama  # or Groq, OpenAI, etc.
# from tools.normalization import normalisation_tool

def main():
    task = input("Enter your analysis task: ")
    data_path = input("Enter your data file path: ")

    # ---- load CSV safely as AnnData & summarize ----
    data = load_csv_as_adata(data_path)
    print("Data shape:", data.shape)
    print("Obs names (first 5):", list(data.obs_names[:5]))
    print("Var names (first 5):", list(data.var_names[:5]))

    data_signature = f"AnnData with {data.n_obs} obs × {data.n_vars} vars"

    # Initialize agents
    llm = Ollama(model='llama3.1', base_url='http://localhost:11430')  # plug Groq/OpenAI later
    memory = GlobalMemory()

    # >>> FIXED: pass data_representation into Planner
    planner = Planner(llm=llm, data_representation=data_signature)

    # >>> If your Planner.plan signature is plan(user_task)
    steps = planner.plan(task)

    print("Planned Steps:", steps)

    for step in steps:
        print(f"Executing step: {step['description']}")
        # if "normalize" in step['description'].lower():
        #     result = normalisation_tool(data_path)
        #     print(result[0].content)

if __name__ == "__main__":
    main()
