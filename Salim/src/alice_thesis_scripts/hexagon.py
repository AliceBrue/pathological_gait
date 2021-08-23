import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.legend_handler import HandlerTuple


def hexagon_plot(txt_file, results_path, param_names, plot_names, place_names, title, param_factors=[1, 1, 1],
                 axes_factor=1):
    """
    Plots hexagonal representation of 3 parameters alteration and pathological gaits
    INPUTS: - txt_file: name of the text file gathering evaluations to plot
            - results_path: path to results folder
            - param_names: list of the 3 parameters of interest
            - plot_names: list of the latex names of the previous parameters for the plot
            - place_names: list of floats to adapt the place of the parameter names on the plot
            - title: plot titles
            - param_factors: optional list of parameter factors to better visualise their alterations
            - axes_factor: otional factor to adapt the scale of the axes
    OUTPUTS: -

    """

    # plot hexagon
    r = 3
    a = np.r_[0, r]
    b = np.r_[-r*np.sqrt(3)/2, r/2]
    c = np.r_[-r*np.sqrt(3)/2, -r/2]
    d = np.r_[0, -r]
    e = np.r_[r*np.sqrt(3)/2, -r/2]
    f = np.r_[r*np.sqrt(3)/2, r/2]
    hexagon = np.c_[a, b, c, d, e , f, a]

    # parameter names at corners
    plot_names = [plot_names[0]+'>0', plot_names[1]+'>0', plot_names[2]+'<0', plot_names[0]+'<0', plot_names[1]+'<0',
                  plot_names[2]+'>0']
    fig, ax = plt.subplots()
    ax.set(xlim=(-r-1.6, r+1.6), ylim=(-r-0.55, r+1.5))
    ax.plot(hexagon[0], hexagon[1], 'k', lw=2)
    ax.plot([a[0], d[0]], [a[1], d[1]], '--k')
    ax.plot([b[0], e[0]], [b[1], e[1]], '--k')
    ax.plot([c[0], f[0]], [c[1], f[1]], '--k')
    ax.text(a[0]+place_names[0], a[1]+0.15, plot_names[0], fontsize=18)
    ax.text(f[0]+place_names[1], f[1]-0.1, plot_names[1], fontsize=18)
    ax.text(e[0]+place_names[2], e[1]-0.3, plot_names[2], fontsize=18)
    ax.text(d[0]+place_names[3], d[1]-0.4, plot_names[3], fontsize=18)
    ax.text(c[0]+place_names[4], c[1]-0.3, plot_names[4], fontsize=18)
    ax.text(b[0]+place_names[5], b[1], plot_names[5], fontsize=18)

    # axes scale
    axes_scale = np.array([-2, -1, 1, 2])
    point_names = [str(i/axes_factor) for i in axes_scale]
    ax.annotate('0', xy=(0, 0), xytext=(0.05, 0.2), fontsize=15)
    ax.plot([0, 0, 0, 0], axes_scale, "_k")
    ax.plot(axes_scale*np.sqrt(3)/2, axes_scale*1/2, "|k")
    ax.plot(-axes_scale*np.sqrt(3)/2, axes_scale*1/2, "|k")
    for i in range(len(axes_scale)):
        ax.annotate(point_names[i], xy=(0.0, axes_scale[i]), xytext=(0.1, axes_scale[i]-0.1), fontsize=15)
        ax.annotate(point_names[i], xy=(axes_scale[i]*np.sqrt(3)/2, axes_scale[i]/2),
                    xytext=(axes_scale[i]*np.sqrt(3)/2 - 0.15, axes_scale[i]/2 + 0.15), fontsize=15)
        ax.annotate(point_names[len(axes_scale)-1-i], xy=(axes_scale[i]*np.sqrt(3)/2, -axes_scale[i]/2),
                    xytext=(axes_scale[i]*np.sqrt(3)/2-0.15, -axes_scale[i]/2+0.15), fontsize=15)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)

    if txt_file == 'example':
        p1, = plt.plot([0], [1], 'Db', markersize=10)
        p2, = plt.plot([1.5*np.sqrt(3)/2], [1.5/2], 'Db', markersize=10)
        p3, = plt.plot([-0.5*np.sqrt(3)/2], [0.5/2], 'Db', markersize=10)
        p4, = plt.plot([0.0], [2.4], 'Dr', markersize=10)
        p5, = plt.plot([1.8*np.sqrt(3)/2], [1.8/2], 'Dr', markersize=10)
        p6, = plt.plot([1*np.sqrt(3)/2], [-1/2], 'D', color='limegreen', markersize=10)
        plt.legend(((p1), p4, p6), ('healthy gait', 'toe gait', 'heel gait'),
                   handler_map={tuple: HandlerTuple(ndivide=None)})

    else:
        file_name = results_path + 'evaluations/' + txt_file + '.txt'
        file = open(file_name, "r")

        n_param = int(file.readline().split('=')[1])
        n_eval = int(file.readline().split('=')[1]) + 1
        var_tab = np.zeros((n_eval, n_param))

        var_names = np.asarray(file.readline().split()[2:])

        n_toe = 0
        n_heel = 0
        n_falling = 0
        line = file.readline()
        for i in range(n_eval):
            if np.asarray(line.split())[1] == 'toe':
                n_toe += 1
            elif np.asarray(line.split())[1] == 'heel':
                n_heel += 1
            elif np.asarray(line.split())[1] == 'falling':
                n_falling += 1
            var_tab[i, :] = np.asarray(line.split())[2:].astype(np.float)
            line = file.readline()

        param_1 = var_tab[:, np.where(var_names == param_names[0])[0][0]]*axes_factor*param_factors[0]
        param_2 = var_tab[:, np.where(var_names == param_names[1])[0][0]]*axes_factor*param_factors[1]
        param_3 = var_tab[:, np.where(var_names == param_names[2])[0][0]]*axes_factor*param_factors[2]

        # plot healthy condition
        p1, = plt.plot([0], param_1[0], 'Db', markersize=10)
        p2, = plt.plot(param_2[0]*np.sqrt(3)/2, param_2[0]*1/2, 'Db', markersize=10)
        p3, = plt.plot(-param_3[0]*np.sqrt(3)/2, param_3[0]*1/2, 'Db', markersize=10)

        # plot toe conditions
        for i in range(n_toe):
            if param_1[i + 1] != param_1[0]:
                p4, = plt.plot([0], param_1[i + 1], 'Dr', markersize=10)
            elif param_2[i + 1] != param_2[0]:
                p5, = plt.plot(param_2[i + 1]*np.sqrt(3)/2, param_2[i + 1]*1/2, 'Dr', markersize=10)
            elif param_3[i + 1] != param_3[0]:
                p6, = plt.plot(-param_3[i + 1]*np.sqrt(3)/2, param_3[i + 1]*1/2, 'Dr', markersize=10)

        # plot heel conditions
        for i in range(n_heel):
            if param_1[i + n_toe + 1] != param_1[0]:
                p7, = plt.plot([0], param_1[i + n_toe + 1], 'D', color='limegreen', markersize=10)
            elif param_2[i + n_toe + 1] != param_2[0]:
                p8, = plt.plot(param_2[i + n_toe + 1]*np.sqrt(3)/2, param_2[i + n_toe + 1]*1/2, 'D', color='limegreen',
                              markersize=10)
            if param_3[i + n_toe + 1] != param_3[0]:
                p9, = plt.plot(-param_3[i + n_toe + 1]*np.sqrt(3)/2, param_3[i + n_toe + 1]*1/2, 'D', color='limegreen',
                               markersize=10)

        # plot falling condition
        for i in range(n_falling):
            if param_1[i + n_eval - n_falling] != param_1[0]:
                p10, = plt.plot([0], param_1[i + n_eval - n_falling], 'Dk', markersize=10)
            elif param_2[i + n_eval - n_falling] != param_2[0]:
                p11, = plt.plot(param_2[i + n_eval - n_falling]*np.sqrt(3)/2, param_2[i + n_eval - n_falling]*1/2, 'Dk',
                                markersize=10)
            elif param_3[i + n_eval - n_falling] != param_3[0]:
                p12, = plt.plot(-param_3[i + n_eval - n_falling]*np.sqrt(3)/2, param_3[i + n_eval - n_falling]*1/2, 'Dk',
                                markersize=10)

        # legend
        if txt_file == 'sp_sw_hexagon' or txt_file == 'cx_TA_hexagon':
            plt.legend(((p1), (p4, p4), (p10)), ('healthy gait', 'toe gait range ', 'falling'),
                       handler_map={tuple: HandlerTuple(ndivide=None)})
        elif txt_file == 'sp_st_hexagon':
            plt.legend(((p1), (p6, p6), (p9, p9), (p10)), ('healthy gait', 'toe gait range', 'heel gait range', 'falling'),
                       handler_map={tuple: HandlerTuple(ndivide=None)})
        elif txt_file == 'cx_SOL_hexagon':
            plt.legend(((p1), (p7, p7), (p10)), ('healthy gait', 'heel gait range', 'falling'),
                       handler_map={tuple: HandlerTuple(ndivide=None)})

    plt.title(title)
    plt.show()
