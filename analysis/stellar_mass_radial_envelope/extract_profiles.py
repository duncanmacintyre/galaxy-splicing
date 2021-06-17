import numpy as np
import matplotlib.pyplot as plt
import pandas

radii = np.loadtxt("../radius_bins.txt")
radii_names = np.array(["{:.1f}kpc".format(n) for n in radii])

def amap(f, *args, dtype=None):
    if dtype == None:
        dtype = args[0].dtype
    return(np.fromiter((f(*i) for i in zip(*(a.flatten() for a in args),)),
                        dtype=dtype, count=args[0].size).reshape(args[0].shape))

# https://stackoverflow.com/a/2566508/13326516
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# def compute_profiles_helper(filepath, interval, fcns=('min', 'median', 'max'), name=None, str_dtype='U32'):

#     # intervals = np.array(interval, dtype=float).reshape((-1, 2))

#     df = pandas.read_fwf(filepath).set_index('Time')

#     interval = amap(lambda x: find_nearest(df.index, x), intervals)

#     if name is None:
#         name = 'from {:.0f} to {:.0f} Myr'.format(interval[0], interval[1])

#     if names is not False:
#         profiles = df.loc[l[0]:(l[1]+0.001),['Total', *radii_names]].aggregate(fcns, axis=0).T
#     else:
#         profiles = pandas.concat([df.loc[l[0]:(l[1]+0.001),['Total', *radii_names]].aggregate(fcns, axis=0).T for l in intervals], axis=1, copy=False)

#     return(profiles, interval, name)


def compute_profiles(filepath, intervals, fcns=('min', 'median', 'max'), names=None, str_dtype='U32'):

    intervals = np.array(intervals, dtype=float).reshape((-1, 2))
    assert(all(intervals[:,0] <= intervals[:,1]))

    if names is None:
        names = True if intervals.shape[0] > 1 else False

    df = pandas.read_fwf(filepath).set_index('Time')

    intervals = amap(lambda x: find_nearest(df.index, x), intervals)
    

    if names is True or names is False:
        interval_names = np.fromiter(('from {:.0f} to {:.0f} Myr'.format(x[0], x[1]) for x in intervals), count=intervals.shape[0], dtype=str_dtype)
    else:
        interval_names = np.array(names)
        assert(interval_names.size==intervals.shape[0])

    if names is not False:
        assert(np.unique(interval_names).size == intervals.shape[0])
        profiles = pandas.concat([df.loc[l[0]:(l[1]+0.001), radii_names].aggregate(fcns, axis=0).T for l in intervals], axis=1, keys=interval_names, copy=False)
        total = pandas.concat([df.loc[l[0]:(l[1]+0.001), 'Total'].aggregate(fcns) for l in intervals], axis=1, keys=interval_names, copy=False)
    else:
        profiles = pandas.concat([df.loc[l[0]:(l[1]+0.001), radii_names].aggregate(fcns, axis=0).T for l in intervals], axis=1, copy=False)
        total = pandas.concat([df.loc[l[0]:(l[1]+0.001), 'Total'].aggregate(fcns) for l in intervals], axis=1, copy=False)

    return(profiles, total, intervals, interval_names)


def compute_profiles_from_files(filepaths, intervals, names=None, **kwargs):

    if names is None:
        names = True if len(filepaths) > 1 else False

    if names is True:
        profiles_tuple, total_tuple, intervals_tuple, interval_names_tuple = zip(*(compute_profiles(f, i, **kwargs) for f, i in zip(filepaths, intervals)))
        interval_names_array = np.concatenate(interval_names_tuple)
        assert(np.unique(interval_names_array).size == len(intervals_tuple)) # ensure each name is unique
        return(
            pandas.concat(profiles_tuple, axis=1, keys=interval_names_array, copy=False), # profiles
            pandas.concat([x.rename({'Total':s}, axis=1) for s,x in zip(interval_names_array, total_tuple)], axis=1, copy=False), # total mass
            np.array([np.array(*i) for i in intervals_tuple]), # intervals
            interval_names_array) # interval names
    elif names is False:
        profiles_tuple, total_tuple, intervals_tuple, interval_names_tuple = zip(*(compute_profiles(f, i, **kwargs) for f, i in zip(filepaths, intervals)))
        return(
            pandas.concat(profiles_tuple, axis=1, copy=False), # profiles
            pandas.concat(total_tuple, axis=1, copy=False), # total mass
            np.array([np.array(*i) for i in intervals_tuple]), # intervals
            np.concatenate(interval_names_tuple)) # interval names
    else:
        assert(len(np.unique(names)) == len(names)) # ensure each name is unique
        profiles_tuple, total_tuple, intervals_tuple, _ = zip(*(compute_profiles(f, i, **kwargs) for f, i in zip(filepaths, intervals)))
        return(
            pandas.concat(profiles_tuple, axis=1, keys=names, copy=False), # profiles
            pandas.concat([x.rename({'Total':s}, axis=1) for s,x in zip(names, total_tuple)], axis=1, copy=False), # total mass
            np.array([np.array(*i) for i in intervals_tuple]), # intervals
            np.array(names)) # interval names




# ------ beginning of script part

profiles, total, intervals, names = compute_profiles_from_files(
    ['SPT2349_1e6_gf0.5/star_mass_data_mean.txt', 'SPT2349_1e6_gf0.9/star_mass_data_mean.txt',
     'SPT2349_1e6_real01/star_mass_data_mean.txt', 'SPT2349_1e6_real02/star_mass_data_mean.txt',
     'SPT2349_1e6_real03/star_mass_data_mean.txt', 'SPT2349_1e6_real04/star_mass_data_mean.txt'],
    [[700, 800], [700, 800], [700, 800], [700, 800], [700, 800], [700, 800]],
    names=['fgas=0.5', 'fgas=0.9', 'fgas=0.7, real01', 'fgas=0.7, real02', 'fgas=0.7, real03', 'fgas=0.7, real04'])

# column names for fgas=0.7
names7 = ['fgas=0.7, real01', 'fgas=0.7, real02', 'fgas=0.7, real03', 'fgas=0.7, real04']

# mean profile for fgas=0.7
mean = profiles.loc[:, [[s, 'median'] for s in names7]].mean(axis=1)

# print everything out
for x in (profiles, total, intervals, names): print(x)
print(mean)

# --- save the mean profile for fgas=0.7

i = intervals[names == names7[0], :].flatten()
assert((intervals[np.isin(names, names7),0] == i[0]).all() and (intervals[np.isin(names, names7),1] == i[1]).all()) # check that all intervals were the same

export_df = pandas.DataFrame(
    pandas.concat((pandas.Series({'Total':total.loc['median',names7].mean()}), mean), axis=0),
    columns=('from {:.0f} to {:.0f} Myr'.format(*i),))
# save to file
export_df.to_string('radial_profiles_{}.txt'.format('fgas=0.7_mean'), float_format=lambda x: '{:.5g}'.format(x))

# --- plot the profiles

# colours = ['seagreen', 'royalblue', 'm', 'firebrick', 'darkolivegreen', 'saddlebrown']
colours_bak = ['springgreen', 'deepskyblue', 'fuchsia', 'darkslateblue', 'orange', 'peru']

fig, ax = plt.subplots(figsize=(8,5), constrained_layout=True)

for s, c in zip(names, colours_bak):
    ax.fill_between(radii, profiles[s, 'min'], profiles[s, 'max'], label=s, color=c)
for s, in zip(names):
    ax.plot(radii, profiles[s, 'median'], color='k', linewidth=0.5)
ax.plot(radii, mean, linestyle='--', label='fgas=0.7, mean')

ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.6))
ax.set_xlabel('Radius [kpc]')
ax.set_ylabel('Mass contained [M☉]')
ax.set_title('Radial distribution from {:.0f} to {:.0f} Myr'.format(*i))

plt.savefig('radial-profiles.pdf')

plt.show()
