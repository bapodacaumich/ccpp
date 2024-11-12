from station import visualize_station
from utils import plot_path_direct, save_animation, plot_cw_opt_path, set_aspect_equal_3d, plot_packaged_path, plot_viewpoints
from matplotlib import pyplot as plt
import os
import numpy as np


def plotpath(vgd, local, condition='ocp'):
    local_txt = '_local' if local else '_global'
    print( 'Plotting path for', vgd + local_txt, ' condition:', condition)
    # ax = visualize_station(coverage_file= vgd + '_global_coverage.csv', vx=False)
    ax = visualize_station(coverage_file=None, convex=False, vx=False)
    # ax = visualize_station(coverage_file=None, convex=True, vx=False)
    # ax = visualize_station(coverage_file=None, convex=False, vx=False, original=True, ax=ax)
    # ax = visualize_station(coverage_file=condition + '/' + vgd+local_txt+'_cov.csv', convex=False, vx=False, final=True)

    if condition[0] == 'i':
        soln_folder = 'cw_opt_packaged_' + condition.split('_')[1]
    else: # condition is ocp
        soln_folder = 'pf_final'
    ax = plot_packaged_path(soln_folder, vgd + local_txt + '.csv', ax=ax, plot_dir=True)
    ax = plot_viewpoints(ax, vgd, local, connect=True)
    # ax = plot_path_direct('unfiltered_viewpoints', 'unfiltered_viewpoints_' + vgd + '.csv', ax=ax)
    # ax = plot_path_direct('coverage_viewpoint_sets', 'coverage_' + vgd + '_local_vp_set.csv', ax=ax, ordered=False, local=True)
    # ax.legend()
    # ax = plot_path_direct('ordered_viewpoints', vgd + local_txt + '.csv', ax=ax, ordered=True, local=local)
    # ax = plot_cw_opt_path(ax, vgd + local_txt)
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax = plot_path_direct('ordered_viewpoints', vgd + '.csv', ax=ax)
    # plt.savefig('4m_unfiltered.png', dpi=300)
    ax.view_init(elev=30, azim=30)
    ax.dist = 5
    ax.set_axis_off()

    set_aspect_equal_3d(ax)
    # save_animation(ax, vgd + local_txt + '_ocp')


    # dir = os.path.join(os.getcwd(), 'figures', condition)
    # if not os.path.exists(dir): os.makedirs(dir)
    # plt.savefig(os.path.join(os.getcwd(), 'figures', condition, vgd + local_txt + '.png'), dpi=600)
    # plt.savefig(os.path.join(os.getcwd(), 'figures', 'ordered_vps', vgd + local_txt + '.png'), dpi=600)
    # plt.tight_layout()
    plt.show()
    # plt.close('all')

    # debug coverage
    # vp0 = np.array([0.574084,6.236693,6.799296, 0.489475,0.033094,-0.871389])
    # ax.scatter(vp0[0], vp0[1], vp0[2], c='k', alpha=1.0)
    # ax.quiver(vp0[0], vp0[1], vp0[2], vp0[3], vp0[4], vp0[5], length=1)

    # ax.scatter(2.029000,2.821000,2.091000, c='r', alpha=1.0)
    # ax.scatter(3.529000,8.821000,2.5910001, c='r', alpha=1.0)
    # start = [2.029000,2.821000,2.091000]
    # end = [1.029000,2.821000,2.091000]
    # # end = [6.183021,2.812072,0.445231]
    # # end = [3.022275,7.297307,5.561572]
    # ray = np.array([start, end])
    # ax.plot(ray[:,0], ray[:,1], ray[:,2], 'k-')
    # ints = np.array([[1.878773,2.821000,2.091000],[1.760119,2.821000,2.091000],[1.992262,2.821000,2.091000],[1.966639,2.821000,2.091000],[1.604648,2.821000,2.091000],[1.575547,2.821000,2.091000],[1.516622,2.821000,2.091000],[1.900703,2.821000,2.091000],[1.956441,2.821000,2.091000],[1.516622,2.821000,2.091000],[1.654887,2.821000,2.091000],[1.844243,2.821000,2.091000],[1.760136,2.821000,2.091000],[1.604648,2.821000,2.091000],[1.579064,2.821000,2.091000],[1.527835,2.821000,2.091000],[1.674147,2.821000,2.091000],[1.867668,2.821000,2.091000],[1.848995,2.821000,2.091000],[1.822615,2.821000,2.091000],[1.809751,2.821000,2.091000],[1.790757,2.821000,2.091000],[1.790757,2.821000,2.091000],[1.790757,2.821000,2.091000],[1.803961,2.821000,2.091000],[1.822615,2.821000,2.091000],[1.842012,2.821000,2.091000],[1.867668,2.821000,2.091000],[1.881771,2.821000,2.091000],[1.857782,2.821000,2.091000],[1.791786,2.821000,2.091000],[1.814568,2.821000,2.091000],[1.845701,2.821000,2.091000],[1.828557,2.821000,2.091000],[1.803688,2.821000,2.091000],[1.789968,2.821000,2.091000],[1.769116,2.821000,2.091000],[1.826814,2.821000,2.091000],[1.899461,2.821000,2.091000],[1.855649,2.821000,2.091000],[1.790980,2.821000,2.091000],[1.813662,2.821000,2.091000],[1.844705,2.821000,2.091000],[1.844648,2.821000,2.091000],[1.844571,2.821000,2.091000],[1.821377,2.821000,2.091000],[1.786954,2.821000,2.091000],[1.836788,2.821000,2.091000],[1.899526,2.821000,2.091000],[1.899640,2.821000,2.091000],[1.899783,2.821000,2.091000],[1.872169,2.821000,2.091000],[1.833772,2.821000,2.091000],[1.827135,2.821000,2.091000],[1.817693,2.821000,2.091000],[1.812227,2.821000,2.091000],[1.804305,2.821000,2.091000],[1.863102,2.821000,2.091000],[1.955132,2.821000,2.091000],[1.843755,2.821000,2.091000],[1.815306,2.821000,2.091000],[1.751132,2.821000,2.091000],[1.926331,2.821000,2.091000],[1.669160,2.821000,2.091000],[1.654865,2.821000,2.091000],[1.654865,2.821000,2.091000],[1.660121,2.821000,2.091000],[1.744171,2.821000,2.091000],[1.846554,2.821000,2.091000],[2.010139,2.821000,2.091000],[1.992514,2.821000,2.091000],[1.766740,2.821000,2.091000],[1.758250,2.821000,2.091000],[1.896451,2.821000,2.091000],[1.981495,2.821000,2.091000],[1.839370,2.821000,2.091000],[1.710612,2.821000,2.091000],[1.710612,2.821000,2.091000],[1.710612,2.821000,2.091000],[1.864654,2.821000,2.091000],[1.946477,2.821000,2.091000],[1.844251,2.821000,2.091000],[1.877458,2.821000,2.091000],[1.918035,2.821000,2.091000],[1.969176,2.821000,2.091000],[1.991713,2.821000,2.091000],[1.958134,2.821000,2.091000],[1.883120,2.821000,2.091000],[1.769116,2.821000,2.091000],[1.826974,2.821000,2.091000],[1.899783,2.821000,2.091000],[1.927028,2.821000,2.091000],[1.945328,2.821000,2.091000],[1.391594,2.821000,2.091000],[1.390864,2.821000,2.091000],[1.391594,2.821000,2.091000],[1.390864,2.821000,2.091000],[1.391594,2.821000,2.091000],[1.390864,2.821000,2.091000],[1.391594,2.821000,2.091000],[1.390864,2.821000,2.091000],[1.391594,2.821000,2.091000],[1.390864,2.821000,2.091000],[1.392167,2.821000,2.091000],[1.444621,2.821000,2.091000],[1.392167,2.821000,2.091000],[1.392167,2.821000,2.091000],[1.392167,2.821000,2.091000],[1.547580,2.821000,2.091000]])
    # ax.plot(ints[:,0], ints[:,1], ints[:,2], 'rx')
    # print('xlim=', ax.get_xlim())
    # print('ylim=', ax.get_ylim())
    # print('zlim=', ax.get_zlim())
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # plt.show()

if __name__ == "__main__":
    condition = 'ocp'
    # vgds = ['2m', '4m', '8m', '16m']
    vgds = ['4m']
    locals = [True, False]
    for vgd in vgds:
        for local in locals:
            plotpath(vgd, local, condition=condition)
    # plotpath()