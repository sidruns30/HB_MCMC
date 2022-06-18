import corner
import numpy as np
import sys, os, glob
import matplotlib.pyplot as plt


"""
Functions to compute the radii and temperature of the stars
"""
def temperature_model(logM):
    # in solar masses
    M_nodes = [0.1, 0.26, 0.47, 0.59, 0.69, 0.87,
              0.98, 1.085, 1.4, 1.65, 2.0, 2.5, 
              3.0, 4.4, 15., 40.]
    # in log10K
    T_nodes = [3.491, 3.531, 3.547, 3.584, 3.644, 3.712,
              3.745, 3.774, 3.823, 3.863, 3.913, 3.991,
              4.057, 4.182, 4.477, 4.623]
    
    m = np.power(10., logM)
    
    if m <= M_nodes[0]: T = T_nodes[0]
    elif m >= M_nodes[15]: T = T_nodes[15]
    else:
        for j, t in enumerate(T_nodes):
            if (m < M_nodes[j]):
                T = T_nodes[j-1] + (m - M_nodes[j-1]) * (T_nodes[j] - T_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1])
                break
    return T

def radius_model(logM):
    # in solar masses
    M_nodes = [0.07, 0.2, 0.356, 0.655, 0.784, 0.787, 1.377,
               4.4, 15., 40.]
    # in log10 Rsun
    R_nodes = [-0.953, -0.627, -0.423, -0.154, -0.082, -0.087,
            0.295, 0.477, 0.792, 1.041]
    
    m = np.power(10., logM)
    
    if m <= M_nodes[0]: R = R_nodes[0]
    elif m >= M_nodes[9]: R = R_nodes[9]
    else:
        for j, r in enumerate(R_nodes):
            if (m < M_nodes[j]):
                R = R_nodes[j-1] + (m - M_nodes[j-1]) * (R_nodes[j] - R_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1])
                break
    return R


def envelope_temp(logM):
    return 0.0224

def envelope_radius(logM):
    m = np.power(10., logM)
    n = 4.22
    slope = 15.68
    floor=0.01
    corner=1.055
    ceil=0.17
    boundary = 1/(1/ceil+1/(slope*np.power((np.power(m, n) + 
                np.power(corner, n)), (1/n)) - (slope*corner-floor)))
    return boundary


def temperature(logM, alphaT):
    T = np.power(10., temperature_model(logM) + alphaT * envelope_temp(logM))
    return T

def radius(logM, rr):
    R = np.power(10., radius_model(logM) + rr * envelope_radius(logM))
    return R


temperature_model = np.vectorize(temperature_model)
radius_model = np.vectorize(radius_model)

envelope_temp = np.vectorize(envelope_temp)
envelope_radius = np.vectorize(envelope_radius)

temperature = np.vectorize(temperature)
radius = np.vectorize(radius)

def generate_relevant_samples(samples):
    # input samples is the nchain x 24 parameter array
    # we first remvove the period and the likelihood from it
    # relevant quantities and their indices
    # log M1: 2
    # log M2: 3
    # e:      5  
    # inc:    6
    # omega:  7
    # rr1:    9
    # rr2:    10
    # teff1:  19
    # teff2:  20

    R1 = radius(samples[:, 2], samples[:, 9])
    R2 = radius(samples[:, 3], samples[:, 10])
    T1 = temperature(samples[:, 2], samples[:, 19])
    T2 = temperature(samples[:, 3], samples[:, 20])

    samples = np.hstack((np.power(10., samples[:,2:4]), samples[:,5].reshape(-1, 1), 
                180*samples[:,6].reshape(-1, 1)/np.pi, 180*samples[:,7].reshape(-1, 1)/np.pi, 
                R1.reshape(-1, 1), R2.reshape(-1, 1),
                T1.reshape(-1,1), T2.reshape(-1, 1), np.power(10., samples[:, 2] - samples[:, 3]).reshape(-1, 1)))
    labels = ["M1", "M2", "e", "inc", "omega", "R1", "R2", "T1", "T2", "M1/M2"]
    return samples, labels

    
"""
Load the corresponding chain and generate relevant corner plots
args:
    lc_id: TIC of the source
    run_id: the run # of the chain name (default to 0)
    skip_samples: number of samples to disregard during the burn in
    gmag: load chain that includes the gmag data
    color: load chain that includes the color data
    omp: load chains from OMP runs
    orbital: only plot the orbital parameters
returns:
    samples: array with chain samples
    labels: the respective labels of the chain data
"""

def load_chain_dat(lc_id, run_id=0, skip_samples=1000, gmag=False, color=False, omp=True, orbital=False, **kwargs):
    chain_loc = "/scratch/ssolanski/HB_MCMC/data/chains/"
    chain_loc += "chain.%s" % lc_id
    if gmag: chain_loc += "_gmag"
    if color: chain_loc += "_color"
    if omp: chain_loc += "_OMP"
    chain_loc += "_%d.dat" % run_id

    
    print("Loading chain %s" % chain_loc)
    chain = np.genfromtxt(chain_loc, skip_header=skip_samples, skip_footer=1)
    print("Number of samples: %d" % chain.shape[0])
    print(np.average(chain[:, 1]))
    if orbital:
        samples, labels = generate_relevant_samples(chain)
        return samples, labels

    # remove the period from the chain
    samples = np.hstack((chain[:, 2:4], chain[:, 5:]))
    labels = ["M1", "M2", "e", "inc", "omega", "T0", "rr1", "rr2", "mu 1", 
              "mu 2", "tau 1", "tau 2", "ref 1", "ref 2", "alpha_b 1", "alpha_b_2", 
              "TT1", "TT2", "blending", "flux_tune"]

    return samples, labels

"""
Load and return an external chain; does not return any labels on the chain
"""
def load_custom_chain(orbital, chain_path):
    data = np.genfromtxt(chain_path, skip_header=10)
    print(data.shape)
    # shape for John's samples
    if data.shape[1] == 23:
        data = data[:, 1:-1]
        print(data[-1, :])
    elif data.shape[1] == 27:
        data = data[:, 4:-2]
        print(data[-1, :])
    else:
        print("Incorrect data shape: ", data.shape)

    if orbital:
        data, labels = generate_relevant_samples(data)
        return data
    else:
        return data[:, 1:]


"""
Collect all the chain data in one list
args: lc_id: TIC ID of the source
      run_ids: all the run ids required for the 
      **load_chain_dat_kwargs: other args to specify
      what chains to load
returns: list of all the chain data
"""

def load_all_data(lc_id, run_ids=[0], **kwargs):
    plot_range = []
    samples_list = []
    for index in run_ids:
        #for gmag in [False, True]:       
        #    kwargs["gmag"] = gmag
        samples, labels = load_chain_dat(lc_id=lc_id, run_id=index, **kwargs)
        samples_list.append(samples)
    if "chain_file_locs" in kwargs:
        for loc in kwargs["chain_file_locs"]:
            samples_list.append(load_custom_chain(chain_path=loc,orbital=kwargs["orbital"]))
    return samples_list, labels

"""
Set the corner plot args and limits based on the sample data
args:  samples_list: list containing all the chain data in the lc
       fig: matplotlib figure on which to create the corner plot, default is none
returns: dict containing the corner kwargs
"""

def set_corner_kwargs(samples_list, fig=None):
    CORNER_KWARGS = dict(
        smooth=None,
        label_kwargs=dict(fontsize=22, weight="light"),
        title_kwargs=dict(fontsize=22, weight="light"),
        title_fmt = ".4f",
        quantiles=[0.16, 0.84],
        levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
        plot_density=False,
        plot_datapoints=True,
        fill_contours=True,
        show_titles=True,
        max_n_ticks=3,
        use_math_text = True
    )

    samples = samples_list[0]

    plot_range = []
    for dim in range(samples.shape[-1]):
        plot_range.append(
            [
                min([min(samples_list[i].T[dim]) for i in range(len(samples_list))]),
                max([max(samples_list[i].T[dim]) for i in range(len(samples_list))]),
            ]
        )
    
    CORNER_KWARGS.update(range=plot_range)
    if fig is not None: CORNER_KWARGS.update(fig=fig)
    return CORNER_KWARGS

"""
Make the corner plot from the input samples and labels
args: samples: the chain samples
      labels: the respective sample labels
      **CORNER_KWARGS: additional corner plot arguments
returns: fig - the mpl figure on which more corner plots can be drawn
"""
def plot_corner(samples, labels, CORNER_KWARGS):
    corner.corner(samples, labels=labels, **CORNER_KWARGS)
    plt.gca().tick_params(axis='both', labelsize=15)
    return plt.gcf()


def main():
    lc_id = sys.argv[1]
    run_ids = [0, 1]
    gmag = False
    color = False
    omp = True
    orbital = False
    skip_samples = 1000

    if gmag:
        chain_file_locs = glob.glob("/scratch/ssolanski/HB_MCMC/data/chains/237753600_cmag*")
    else:
        chain_file_locs = glob.glob("/scratch/ssolanski/HB_MCMC/data/chains/237753600_std**")

    samples_list, labels = load_all_data(lc_id, run_ids, 
                            gmag=gmag, color=color, omp=omp, orbital=orbital,
                            skip_samples=skip_samples, chain_file_locs=chain_file_locs)

    colors_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
    legends_list = ["TIC_%s_OMP_%s_GMAG_%s_COLOR_%s_ID_%d" % (str(lc_id), str(omp), str(gmag), str(color), run_id) for run_id in run_ids]
    for loc in chain_file_locs: legends_list.append(loc)
    CORNER_KWARGS = set_corner_kwargs(samples_list)

    for i, samples in enumerate(samples_list):
        fig = plot_corner(samples, labels, CORNER_KWARGS)
        CORNER_KWARGS.update(fig=fig)
        CORNER_KWARGS.update(color=colors_list[i])
        CORNER_KWARGS.update(legend=legends_list[i])
        fig.gca().legend()

    plot_name = "../extra/corner_%s" % lc_id
    if gmag: plot_name += "_gmag"
    if color: plot_name += "_color"
    if omp: plot_name += "_omp"
    if orbital: plot_name += "_orbital"
    plot_name += "_all.png"

    plt.savefig(plot_name, dpi=100)

if __name__ == "__main__":
    main()
