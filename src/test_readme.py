



from biograpy import Panel, tracks, features
panel=Panel(fig_width=1000, )#initialize a drawer
test_track = tracks.BaseTrack(features.Simple(name='feat1',start=100,end=756,),
            features.Simple(name='feat2',start=300,end=1056,),
            features.Simple(name='feat3',start=600,end=1356,),
            features.Simple(name='feat4',start=800,end=1356,),
            features.Simple(name='feat5',start= 1357,end=1806,),
            name = 'test')
panel.append(test_track)

panel.save('test.png')
