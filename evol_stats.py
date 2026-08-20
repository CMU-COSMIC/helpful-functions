'''
Writen by Hannah Skobe
'''

###################################################################################################
# natal kicks
###################################################################################################

def split_kicks(kick_df, lower, upper):
    '''
    Separates the kick dataframe into three dataframes based on the chosen lower and upper kick 
    magnitudes for the primary star.

    Parameters
    ----------

    kick_df : pandas.Dataframe
        Returned kick_info dataframe from COSMIC or BackPop
    
    lower : int or float
        A lower kick magnitude (km/s) to set low kick boundary
    
   upper : int or float
        An upper kick magnitude (km/s) to set high kick boundary   

    Returns
    -------

    low : pandas.Datafram
        A dataframe of kick_info for the primary star with an upper bound of the lower input
    
    med : pandas.Datafram
        A dataframe of kick_info for the primary star with a lower bound of the lower input and
        upper bound of the upper input
    
    high : pandas.Datafram
        A dataframe of kick_info for the primary star with a lower bound of the upper input

    '''

    low = kick_df.loc[(kick_df.natal_kick < lower) & (kick_df.star == 1.0)]
    med = kick_df.loc[(kick_df.natal_kick > upper) & (kick_df.natal_kick < upper) & (kick_df.star == 1.0)]
    high = kick_df.loc[(kick_df.natal_kick > upper) & (kick_df.star == 1.0)]
    
    return low, med, high

###################################################################################################
# evolutionary steps
###################################################################################################

def evol_steps(df):
    '''
    Separates a binary population dataframe into key evolutionary steps, namely: the start/end of 
    the first Common Envolope, the start/end of the first Roche-Lobe Overflow, and the start/end 
    of the first supernova of the primary.
    
    Parameters
    ----------

    df : pandas.DataFrame
        A dataframe (e.g., a bpp from COSMIC or BackPop) of binaries
    
    Returns
    -------

    CE1 : pandas.DataFrame
        A dataframe of the first instance of the start of Common Envelope (i.e., evol_type = 7)

    post_CE1 : pandas.DataFrame
        A dataframe of the first instance of the end of Common Envelope (i.e., evol_type = 8)

    RLO1 : pandas.DataFrame
        A dataframe of the first instance of the start of Roche-Lobe Overflow (i.e., evol_type = 3)

    post_RLO1 : pandas.DataFrame
        A dataframe of the first instance of the end of Roche-Lobe Overflow (i.e., evol_type = 4)

    RLO2 : pandas.DataFrame
            A dataframe of the second instance of the start of Roche-Lobe Overflow (i.e., evol_type = 3)

    SN1 : pandas.DataFrame
        A dataframe of the start of the first supernova of the primary (i.e., evol_type = 15)

    post_SN1 : pandas.DataFrame
        A dataframe of the end of the first supernova of the primary (i.e., kstar_1 = 14)
    
    '''

    evol_type_7 = df.loc[(df.evol_type == 7)].groupby('bin_num', as_index=False).first()
    post_CE1 = df.loc[(df.evol_type == 8)].groupby('bin_num', as_index=False).first()
    CE1 = evol_type_7.loc[(evol_type_7.bin_num.isin(post_CE1.bin_num))] # ensure all CE have start and end

    evol_type_3 = df.loc[(df.evol_type == 3)].groupby('bin_num', as_index=False).first()
    post_RLO1 = df.loc[(df.evol_type == 4)].groupby('bin_num', as_index=False).first()
    RLO1 = evol_type_3.loc[(evol_type_3.bin_num.isin(post_RLO1.bin_num))] # ensure all RLO have start and end
    RLO2 = evol_type_3.loc[~(evol_type_3.bin_num.isin(post_RLO1.bin_num))]

    SN1 = df.loc[(df.evol_type == 15)].groupby('bin_num', as_index=False).first()
    post_SN1 = df.loc[(df.kstar_1 == 14)].groupby('bin_num', as_index=False).first()
    
    return CE1, post_CE1, RLO1, post_RLO1, RLO2, SN1, post_SN1


###################################################################################################
# evolutionary stats
###################################################################################################

def evol_stats(df, kick_df):
    '''
    Prints the evolutionary statistics of a binary population dataframe and a kick dataframe:
        - Average initial mass of the primary and secondary
        - Average time spent in the first RLO
        - Number of binaries that transition from SMT to CE before first SN
        - Number of binaries that evolve through CE only before first SN
        - Number of binaries that evolve through SMT only before first SN
        - Number of binaries that evolve through no mass transfer before first SN
        - The primary star type at start of the first RLO
        - Eccentricity before/after the first SN
        - Orbital Period before/after the first SN
        - Mass lost due to the first SN
        - Average natal kick
        - Average mass of the primary and secondary after the first SN
        - Semi-latus rectum after the first SN

    Parameters
    ----------

    df : pandas.DataFrame
        A dataframe (e.g., bpp from COSMIC or BackPop) of binaries
    
    kick_info : 
        Returned kick_info dataframe from COSMIC or BackPop

    '''
    
    CE1, post_CE1, RLO1, post_RLO1, RLO2, SN1, post_SN1 = evol_steps(df)
    
    ########################
    # Initial masses
    ########################

    avg_init_m1 = df.loc[(df.tphys == 0)]['mass_1'].values.mean()
    avg_init_m2 = df.loc[(df.tphys == 0)]['mass_2'].values.mean()

    print(f'Initial m1 mass: {avg_init_m1:.3f}')
    print(f'Initial m2 mass: {avg_init_m2:.3f}')
    print('')

    ########################
    # Mass transfer
    ########################

    t_rlo1 = RLO1.tphys.values
    t_end_rlo1 = post_RLO1.tphys.values
    
    RLO1_time = (t_end_rlo1-t_rlo1).mean()
    
    t_CE1 = CE1.tphys.values
    
    print(f'Avg time in RLO1: {RLO1_time:.3f} Myr')
    print('')
    
    t_rlo_CE = RLO1.loc[(RLO1.bin_num.isin(CE1.bin_num))].tphys.values

    t_transition = t_CE1-t_rlo_CE
    t_transition_frac = len(t_transition[t_transition > 0])/len(RLO1)
    SMT_CE = len(t_transition[t_transition > 0])

    RLO1_frac = len(RLO1.bin_num.unique())/len(df.bin_num.unique())
    RLO1_only = len(RLO1.bin_num.unique())-len(t_transition[t_transition>0])
    
    CE1_frac = len(CE1.bin_num.unique())/len(df.bin_num.unique())
    CE1_only = len(CE1.bin_num.unique())-len(t_transition[t_transition>0])
    
    No_MT_frac = len(RLO2.bin_num.unique())/len(df.bin_num.unique())
    No_MT = len(RLO2.bin_num.unique())
    
    print(f'No. of SMT -> CE: {SMT_CE}, {t_transition_frac*100:.3f}%')
    print(f'No. of   only CE: {CE1_only}, {(CE1_frac-t_transition_frac)*100:.3f}%')
    print(f'No. of  only SMT: {RLO1_only}, {(RLO1_frac-t_transition_frac)*100:.3f}%')
    print(f'No. of     no MT: {No_MT}, {No_MT_frac*100:.3f}%')
    print('')

    ########################
    # kstar_1 at RLO1
    ########################

    for kstar_1_val in RLO1.kstar_1.unique():
        num_kstar_1 = RLO1.loc[(RLO1.kstar_1==kstar_1_val)]
        frac_kstar_1 = len(num_kstar_1)/len(RLO1)

        print(f'No. of kstar_1 = {kstar_1_val} at RLO1: {len(num_kstar_1)}, {frac_kstar_1:.4f}%')
    print('')
    
    ########################
    # ecc
    ########################

    print(f'ecc   at SN: {SN1.ecc.mean():.2f}')
    print(f'ecc post-SN: {post_SN1.ecc.mean():.2f}')
    print('')

    ########################
    # porb
    ########################

    print(f'orbital period   at SN: {SN1.porb.mean():.2f} days')
    print(f'orbital period post-SN: {post_SN1.porb.mean():.2f} days')
    print('')

    ########################
    # SN mass loss
    ########################

    dm1 = SN1.mass_1 - post_SN1.mass_1

    print(f'Mass loss due to SN: {dm1.mean():.2f} M')
    print('')

    ########################
    # natal kick
    ########################

    avg_kick = kick_df['natal_kick'].values.mean()

    print(f'Natal kick: {avg_kick:.2f} km/s')
    print('')

    ########################
    # Final masses
    ########################

    avg_m1 = post_SN1.mass_1.mean()
    avg_m2 = post_SN1.mass_2.mean()

    print(f'Post SN1 m1 mass: {avg_m1:.3f} M')
    print(f'Post SN1 m2 mass: {avg_m2:.3f} M')
    print('')

    ########################
    # semi-latus rectum
    ########################

    slr = post_SN1.sep.mean()*(1 - post_SN1.ecc.mean()**2)

    print(f'a(1-e**2): {slr:.3f}')
    print('')