from cProfile import label
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import argparse


plt.style.use(["science", "muted"])


def plot_distribution_density(data, name, axs):
    kde_factor = 5e-2
    kernel_density = stats.gaussian_kde(data, bw_method=kde_factor)
    x = np.linspace(data.min(), data.max(), 5000)
    plot_density = kernel_density(x)
    axs.plot(x, plot_density, label=name[-5:-4:])
    axs.legend()
    
def plot_file(file_name, axs):
    data = np.loadtxt(file_name)
    axs.plot(data, label=file_name[-5:-4:])
    
    
  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the distribution of a given data set")
    parser.add_argument("-f", "--file", help="The file to read the data from", action="append", nargs="+", dest="files")
    args = parser.parse_args()
    list_files = args.files[0]
    # fig, axs = plt.subplots(1, 1, figsize=(3.5, 3))
    # for file in list_files:
    #     plot_file(file, axs)
    # axs.set_yscale("log")
    
    fig, axs = plt.subplots(1, 1, figsize=(3.5, 3))
    list_size = []
    list_means = []
    list_std = []
    for file in list_files:
        data = np.loadtxt(file)
        list_size.append(int(file[-5:-4:]))
        list_means.append(np.mean(data))
        list_std.append(np.std(data))
        plot_distribution_density(data, file, axs)
    axs.set_xlabel("Determinants")
    axs.set_yscale("log")
    # axs.set_xlim(-8, 8)
    fig.savefig("distribution_dets.pdf")
    
    fig2, axs2 = plt.subplots(2, 1, figsize=(3.5, 3))
    axs2[0].plot(list_size, list_means, "-o")
    axs2[1].plot(list_size, list_std, "-o")
    plt.show()
    plt.close()
    exit(0)

  
