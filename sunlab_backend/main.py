from agents.planner import Planner
from agents.executor import Executor
from agents.evaluator import Evaluator
from agents.memory import GlobalMemory
import scanpy as sc
from tools.normalization import normalisation_tool

def main():
    task = input("Enter your analysis task: ")
    data_path = input("Enter your data file path: ")

    # Summarize data
    data = sc.read_csv(data_path)
    data_signature = str(data)

    # Initialize agents
    llm = None  # You can use Groq/OpenAI here
    memory = GlobalMemory()
    planner = Planner(llm)
    executor = Executor(llm, data_path, memory)
    evaluator = Evaluator(llm)

    steps = planner.plan(task, data_signature)
    print("Planned Steps:", steps)

    for step in steps:
        print(f"Executing step: {step['description']}")
        if "normalize" in step['description'].lower():
            result = normalisation_tool(data_path)
            print(result[0].content)

if __name__ == "__main__":
    main()
