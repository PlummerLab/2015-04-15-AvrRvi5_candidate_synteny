fig = plt.figure(figsize=(6, 9))
gs = gridspec.GridSpec(3, 1, hspace=1, height_ratios=[3, 2, 2])
g1 = compound_axis(
    [[0, 100], [0, 100]],
    [[0, 1], [0, 40], [1, 2]],
    gs[0],
    fig,
    xscaling=300,
    yscaling=[1, 1, 1],
    hspace=0,
    )

g2 = compound_axis(
    [[0, 100], [100, 200], [200, 300]],
    [[0, 1], [0, 40]],
    gs[1],
    fig,
    xscaling=300,
    hspace=0,
    )

g3 = compound_axis(
    [[0, 100], [100, 200], [200, 300]],
    [[0, 1], [0, 40]],
    gs[2],
    fig,
    xscaling=300,
    hspace=0,
    )

for ax in np.concatenate([g1[0], g2[0], g3[0]]):
    if ax is not None:
        ax.xaxis.set_major_formatter(SeqFormatter())

scaffold = 'scaffold_996'

track = FeatureTrack()
feats = [f for f in GENOMES['MNH120'][scaffold].features if f.type == 'CDS'][4:6]

intron = new_shape(Rectangle, width=0.2, y_offset=-0.1, facecolor='black')
exon = new_shape(Rectangle, width=1, y_offset=-0.5, facecolor='black', linewidth=0.)
last_exon = new_shape(Triangle, width=1, y_offset=-0.5, facecolor='black')

CDS = new_shape(
    Feature,
    shape=exon,
    last_shape=last_exon,
    )

track.add_bpfeatures(feats, obj=CDS)

#feature_patches = track.draw(g1[0, 1])

link = new_shape(CrossLink, ax1=g2[1,0], ax2=g2[1, 2], ax1_yrange=(2, 0), ax1_cpoint=-5, ax2_cpoint=-5, facecolor='red', linewidth=0., alpha=0.5)
links = LinkCollection(link, add_to_fig=True)
links.add([[20, 60], [220, 240]])
links.draw()

link = new_shape(CrossLink, ax1=g2[1,0], ax1_yrange=(2, 0), ax1_cpoint=-5, ax2_cpoint=-5, facecolor='green', linewidth=0., alpha=0.5)
links = LinkCollection(link, add_to_fig=True)
links.add([[10, 20], [80, 100]])
links.draw()

link = new_shape(CrossLink, ax1=g2[0,0], ax2=g1[2, 1], ax1_yrange=(-1, 1), ax2_yrange=(3, 0), ax1_cpoint=2, facecolor='blue', linewidth=0., alpha=0.5)
links = LinkCollection(link, add_to_fig=True)
links.add([[60, 10], [80, 60]])
links.draw()


link = new_shape(CrossLink, ax1=g3[0,2], ax2=g1[2, 1], ax1_yrange=(-1, 1), ax2_yrange=(3, 0), facecolor='orange', linewidth=0., alpha=0.5)
links = LinkCollection(link, add_to_fig=True)
links.add([[210, 260], [80, 60]])
links.draw()
