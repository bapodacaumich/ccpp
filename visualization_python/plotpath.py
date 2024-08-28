from station import visualize_station
from utils import plot_path_direct
from matplotlib import pyplot as plt

def plotpath():
    ax = visualize_station(coverage_file='coverage_4m_coverage.csv')
    ax = plot_path_direct('coverage_4m_vp_set.csv', ax=ax)
    plt.show()

if __name__ == "__main__":
    plotpath()