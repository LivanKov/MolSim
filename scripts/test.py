import matplotlib.pyplot as plt
import numpy as np
import platform
import psutil

def plot_benchmarks(labels, results, title="Benchmark Results", xlabel="Options", ylabel="Execution Time (s)", save_path=None):
    x = np.arange(len(labels))
    width = 0.5
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(x, results, width, color='skyblue', edgecolor='black')
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    
    # Label bars with values
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha='center', va='bottom')
    
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Get system specs
    specs = get_system_specs()
    plt.figtext(0.5, -0.1, specs, wrap=True, horizontalalignment='center', fontsize=10)
    
    if save_path:
        plt.savefig(save_path, format=save_path.split('.')[-1], bbox_inches='tight')
    else:
        plt.show()

def get_system_specs():
    cpu_name = platform.processor()
    ram_size = round(psutil.virtual_memory().total / (1024 ** 3), 2)
    os_info = f"{platform.system()} {platform.release()}"
    return f"System: {os_info} | CPU: {cpu_name} | RAM: {ram_size} GB"

# Example usage
labels = ["Option 1", "Option 2", "Option 3"]
results = [1.23, 2.45, 1.89]  # Replace with actual benchmark results
plot_benchmarks(labels, results, save_path="../assets/multithreading_benchmark_results.svg")

# Installation Instructions
# To install the required packages, run the following command:
# pip install matplotlib numpy psutil

