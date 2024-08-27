from station import visualize_station
from utils import plot_path_direct
from matplotlib import pyplot as plt

def plotpath():
    ax = visualize_station()
    ax = plot_path_direct('coverage_2m.csv', ax=ax)
    plt.show()

if __name__ == "__main__":
    plotpath()