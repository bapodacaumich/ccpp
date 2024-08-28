from station import visualize_station
from utils import plot_path_direct
from matplotlib import pyplot as plt

def plotpath():
    ax = visualize_station()
    ax = plot_path_direct('coverage_4m.csv', ax=ax)
    print('xlim=', ax.get_xlim())
    print('ylim=', ax.get_ylim())
    print('zlim=', ax.get_zlim())
    plt.show()

if __name__ == "__main__":
    plotpath()