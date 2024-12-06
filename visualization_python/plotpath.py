from station import visualize_station, station_monotone, station_saturation
from utils import plot_path_direct, save_animation, plot_cw_opt_path, set_aspect_equal_3d, plot_packaged_path, plot_viewpoints, camera_fov_points, get_hex_color_tableau
from matplotlib import pyplot as plt
import os
import numpy as np
import plotly.graph_objects as go

def plotly_add_camera_fov(figure, data, step=5):
    for i in range(0, data.shape[0], step):
        viewdir_pts = camera_fov_points(data[i,:3], data[i,3:6])
        figure.add_trace(go.Scatter3d(
            x=viewdir_pts[:,0], y=viewdir_pts[:,1], z=viewdir_pts[:,2],
            mode='lines',
            line=dict(color='black', width=2),
            showlegend=False
        ))

    if data.shape[0] % step != 1:
        i = data.shape[0] - 1
        viewdir_pts = camera_fov_points(data[i,:3], data[i,3:6])
        figure.add_trace(go.Scatter3d(
            x=viewdir_pts[:,0], y=viewdir_pts[:,1], z=viewdir_pts[:,2],
            mode='lines',
            line=dict(color='black', width=2),
            name='Camera FOV'
        ))

    return figure

def plot_path_and_dir(data, figure, savefile=os.path.join('figures','path'), save=False, show=True):
    """plot path and direction in plotly

    Args:
        data (np.ndarray): N x 6 array of [x, y, z, u, v, w] data
        figure (_type_): Plotly figure
        savefile (str, optional): file string to save html of figure to. Defaults to os.path.join('figures','path').
        save (bool, optional): save to savefile?. Defaults to False.
        show (bool, optional): plot figure in browser. Defaults to True.

    Returns:
        _type_: _description_
    """

    figure.add_trace(go.Scatter3d(
        x=data[:,0], y=data[:,1], z=data[:,2],
        mode='lines', name='Path',
        line=dict(color='black', width=5)
    ))

    figure = plotly_add_camera_fov(figure, data, step=5)

    if save:
        savefile = os.path.join(os.getcwd(), savefile)
        os.makedirs(os.path.dirname(savefile), exist_ok=True)
        figure.write_html(savefile + '.html')

    if show: figure.show()

    return figure

def plot_final_path(vgd, local, condition='ocp'):
    """plot final path in plotly

    Args:
        vgd (str): viewpoint generation distance (i.e. '2m', '4m', '8m', '16m')
        local (bool): local or global viewpoint generation
        condition (str, optional): path generation method (folder within ./final_paths/ directory). Defaults to 'ocp'.
    """

    figure = station_monotone(True, title=condition + ': ' + vgd + (' local' if local else ' global') + ' path', save=False, show=False)
    # figure = station_monotone(False, title=condition + ': ' + vgd + (' local' if local else ' global') + ' path', save=False, show=False, fig=figure)
    data = np.loadtxt(os.path.join('final_paths', condition, vgd + '_local.csv' if local else vgd + '_global.csv'), delimiter=',')

    figure.update_layout(template='simple_white')

    save_file = os.path.join(os.getcwd(), 'figures', 'traj', condition, vgd + ('_local' if local else '_global'))
    os.makdirs(os.path.dirname(save_file), exist_ok=True)
    figure = plot_path_and_dir(data, figure, savefile=save_file, save=True, show=False)

def plot_final_path_min_aoi(vgd, local, condition='ocp'):
    """plot final path with min aoi displayed on station in plotly, then save to figures folder

    Args:
        vgd (_type_): _description_
        local (_type_): _description_
        condition (str, optional): _description_. Defaults to 'ocp'.
    """
    figure = station_saturation(condition, vgd + ('_local' if local else '_global'), 'min', save=False, show=False)

    data = np.loadtxt(os.path.join('final_paths', condition, vgd + '_local.csv' if local else vgd + '_global.csv'), delimiter=',')

    figure.update_layout(template='simple_white')

    save_file = os.path.join(os.getcwd(), 'figures', 'traj_min_aoi', condition, vgd + ('_local' if local else '_global') + '_min_aoi')
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    figure = plot_path_and_dir(data, figure, savefile=save_file, save=True, show=False)

def plot_aoi_condition(condition):
    """plot and save final paths in plotly with min aoi data

    Args:
        condition (str): path generation method (folder within ./final_paths/ directory)
    """

    vgds = ['2m', '4m', '8m', '16m']
    locals = [True, False]
    for vgd in vgds:
        for local in locals:
            plot_final_path_min_aoi(vgd, local, condition=condition)

def plot_condition(condition):
    """plot final paths in plotly

    Args:
        condition (str): path generation method (folder within ./final_paths/ directory)
    """

    vgds = ['2m', '4m', '8m', '16m']
    locals = [True, False]
    for vgd in vgds:
        for local in locals:
            plot_final_path(vgd, local, condition=condition)

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

def plot_all_aoi_traj():
    conditions = [
        'ivt_10',
        'ivt_50',
        'ivt_var',
        'ocp_ko',
        'ocp_ko_slerp',
        'ocp_ma',
        'ocp_ma_slerp',
        'ocp_so'
    ]

    for c in conditions:
        plot_aoi_condition(c)

if __name__ == "__main__":
    plot_all_aoi_traj()
    # plot_condition('ocp_station_oriented')
    # condition = 'ocp'
    # # vgds = ['2m', '4m', '8m', '16m']
    # vgds = ['4m']
    # locals = [True, False]
    # for vgd in vgds:
    #     for local in locals:
    #         plotpath(vgd, local, condition=condition)
    # # plotpath()