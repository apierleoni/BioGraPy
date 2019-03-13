'''
Created on 30/ott/2009

@author: andreapierleoni
'''
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from biograpy.drawer import Panel
from biograpy import features, tracks


class SeqRecordDrawer(object):
    '''
    Takes in input a seqrecord and draws all the included features along a 
    reference object
    
    
        ======================= ================================================
        Key                      Description
        ======================= ================================================
        draw_type               ``'simple'`` | ``'gene structure'``, default is
                                ``'simple'``. if `simple` draw is selected all 
                                the features will be drawn. 
                                if `gene structure` is selected only features 
                                reporting to a given gene in 
                                ``feature.qualifeirs['gene']`` will be drawn. 
                                Gene can be set using the `gene_to_draw` 
                                function.
        gene_to_draw            default is ``None`` if a string is specified 
                                only features  carrying that string in the
                                ``feature.qualifeirs['gene']`` qualifier will
                                be drawn.  if ``None``, the first feature of 
                                type ``'gene'``is selected as the gene to draw.
        feature_class_qualifier assign the qualifier key to use to set a 
                                feature_class
        feature_name_qualifier  assign the qualifier key to use to set a 
                                feature_name                        
                                                 
        ======================= ================================================


    all supplied kwargs are passed to the Panel object that is available at
    self.panel

    '''


    def __init__(self,
                 seqrecord,
                 draw_type='simple', 
                 gene_to_draw=None,
                 feature_class_qualifier='ft_description', 
                 feature_name_qualifier='hit_name',
                 **kwargs):

        self.seqrecord=seqrecord
        self.panel = Panel(**kwargs)
        self.feature_class_qualifier=feature_class_qualifier
        self.feature_name_qualifier=feature_name_qualifier

        if draw_type=='simple':
            ref_obj_track = tracks.BaseTrack(features.Simple(1,len(self.seqrecord),height= 1.5, name=seqrecord.name,color_by_cm = True, track_lines = 3,  alpha =1,))
            self.panel.add_track(ref_obj_track)
            self.draw_features()
        elif draw_type=='gene structure':
            self.draw_features_ordered_by_gene_struct(gene_code=gene_to_draw)

    def save(self,output,**kwargs):
        self.panel.save(output,**kwargs)
        self.panel.close()

    def imagemap(self, **kwargs):
        return self.panel.imagemap(**kwargs)

    def boxes(self):
        return self.panel.boxes()

    def draw_features(self):
        tracks2draw = dict()
        feat_collector={}#store complex feature grouping. drawing for witch it is called later
        for feature in self.seqrecord.features:
            if feature.type not in tracks2draw:
                tracks2draw[feature.type] = tracks.BaseTrack(name = feature.type)
            feature_name = ''
            if feature.id != '<unknown id>':
                feature_name=feature.id
            if self.feature_name_qualifier in feature.qualifiers:
                if isinstance(feature.qualifiers[self.feature_name_qualifier], list):
                    feature_name=feature.qualifiers[self.feature_name_qualifier][0]
                elif isinstance(feature.qualifiers[self.feature_name_qualifier], str):
                    feature_name=feature.qualifiers[self.feature_name_qualifier]

            '''draw segmented features '''
            drawn=False
            if feature.sub_features:
                for sub_feature in feature.sub_features:#detect 'join subfeature and calls the appropriate feature
                    if sub_feature.location_operator=='join':
                        grfeat = features.SegmentedSeqFeature(feature,name=feature_name,)
                        tracks2draw[feature.type].add_feature(grfeat)
                        drawn=True
                        break

            else:
                if not drawn:
                    '''Add here feature specific visualizations '''
                    if feature.type=='gene':#'gene' specific visualization
                        grfeat = features.GeneSeqFeature(feature, name=feature_name,)
                        tracks2draw[feature.type].add_feature(grfeat)
                    elif feature.location.start.position == feature.location.end.position:#single point feature
                        grfeat=features.SinglePositionFeature(feature, name=feature_name,)
                        tracks2draw[feature.type].add_feature(grfeat)
                    elif 'transmembrane' in feature.type.lower(): #add TM handling, cannot handle multiple TM annotations
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'transmembrane' in feat_collector:
                            feat_collector['transmembrane']=dict(TM = [], cyto = [], non_cyto = [], )
                        feat_collector['transmembrane']['TM'].append(feature)
                        if 'name' not in feat_collector['transmembrane']:
                            feat_collector['transmembrane']['name'] = feature_name
                    elif 'topological domain' in feature.type.lower(): #handles topological domains from uniprot
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'transmembrane' in feat_collector:
                            feat_collector['transmembrane']=dict(TM = [], cyto = [], non_cyto = [])
                        if 'description' in feature.qualifiers:
                            cyto=False
                            if isinstance(feature.qualifiers['description'],list):
                                for d in feature.qualifiers['description']:
                                    if d.lower() == 'cytoplasmic':
                                        cyto= True
                            else:
                                if feature.qualifiers['description'].lower() == 'cytoplasmic':
                                     cyto= True
                            if cyto:
                                feat_collector['transmembrane']['cyto'].append(feature)
                            else:
                                feat_collector['transmembrane']['non_cyto'].append(feature)
                    elif 'cytoplasm' in feature.type.lower():#risky!!
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'transmembrane' in feat_collector:
                            feat_collector['transmembrane']=dict(TM = [], cyto = [], non_cyto = [])
                        if ('not' in feature.type.lower()) or ('non' in feature.type.lower()):
                            feat_collector['transmembrane']['non_cyto'].append(feature)
                        else:
                            feat_collector['transmembrane']['cyto'].append(feature)

                    elif feature.type=='beta_strand':#add secondary structure handling
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'SecStructFeature' in feat_collector:
                            feat_collector['secondary structure']=dict(betas = [], alphah = [], coil = [])
                        feat_collector['secondary structure']['betas'].append(feature)
                    elif feature.type=='alpha_helix':#add secondary structure handling
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'secondary structure' in feat_collector:
                            feat_collector['secondary structure']=dict(betas = [], alphah = [], coil = [])
                        feat_collector['secondary structure']['alphah'].append(feature)
                    elif feature.type=='peptide_coil':#add secondary structure handling
                        if feature.type in tracks2draw:
                            del tracks2draw[feature.type]
                        if not 'secondary structure' in feat_collector:
                            feat_collector['secondary structure']=dict(betas = [], alphah = [], coil = [])
                        feat_collector['secondary structure']['coil'].append(feature)
                    else:#generic feature
                        if 'score' in feature.qualifiers:#color by score
                            try:
                                color_by_score = False
                                if isinstance(feature.qualifiers['score'],list):
                                    score = float(feature.qualifiers['score'][0])
                                    if 0 <= score <= 1:
                                        feat_score = score
                                        color_by_score = True
                                else:
                                    score = float(feature.qualifiers['score'])
                                    if 0 <= score <= 1:
                                        feat_score = score
                                        color_by_score = True
                                if color_by_score:
                                    grfeat=features.GenericSeqFeature(feature,
                                                                      name=feature_name, 
                                                                      score = feat_score, 
                                                                      use_score_for_color = True, 
                                                                      cm = 'RdYlGn')
                                else:
                                    raise ValueError('Score value cannot be handled, should be a float between 0 and 1. Actual value: %f'%score)
                            except: #draw normal feature
                                grfeat=features.GenericSeqFeature(feature, name=feature_name)
                        else:
                            grfeat=features.GenericSeqFeature(feature, name=feature_name)
                        tracks2draw[feature.type].add_feature(grfeat)

        '''draw complex features '''
        for feat_group in feat_collector:
            if feat_group not in tracks2draw:
                tracks2draw[feat_group] = tracks.BaseTrack(name = feat_group)
            if feat_group == 'transmembrane':
                grfeat=features.TMFeature(fill_color = 'k',
                                         cyto_color = 'r',
                                         #cyto_label = 'in',
                                         non_cyto_color = 'g',
                                         non_cyto_linestyle = '-',
                                         #non_cyto_label = 'out',
                                         TM_ec = 'orange',
                                         TM_fc = 'orange',
                                         TM_label = '',
                                         **feat_collector[feat_group])
            elif feat_group == 'secondary structure':
                grfeat=features.SecStructFeature(**feat_collector[feat_group])
            tracks2draw[feat_group].add_feature(grfeat)
        '''add tracks to panel'''
        for key in sorted(tracks2draw.keys()):
            self.panel.add_track(tracks2draw[key])
          
        if self.seqrecord.letter_annotations:
            per_letter_feats = []
            ymax = ymin = 0
            for quality in self.seqrecord.letter_annotations:
                if isinstance(self.seqrecord.letter_annotations[quality][0], (int,float)):
                    '''draw numeric per letter annotation as line plot feature '''
                    feat = features.PlotFeature(self.seqrecord.letter_annotations[quality], range(1, len(self.seqrecord.letter_annotations[quality])+1), label = str(quality),) 
                    per_letter_feats.append(feat)
                    if min(self.seqrecord.letter_annotations[quality])<= ymin:
                        ymin = min(self.seqrecord.letter_annotations[quality])
                    if max(self.seqrecord.letter_annotations[quality])>= ymax:
                        ymax = max(self.seqrecord.letter_annotations[quality])
            
            per_letter_track = tracks.PlotTrack(ymin = ymin, ymax = ymax)
            for feat in per_letter_feats:
                per_letter_track.append(feat)
            self.panel.add_track(per_letter_track)        
            
            



    def draw_features_ordered_by_gene_struct(self,gene_code=None):

        '''
        Draws 'gene structure' type.
        The input SeqRecord should include a gene, mRNAs, CDSs and other allowed
        features. they will be drawn in a defined order  after fusing 
        corresponding mRNA and CDS. if exons are passed in the seqrecord they 
        are mapped to the gene patch then corresponding mRNA and CDS are coupled
        and drawn then misc_RNA are drawed in the end the other feature present
        in the seqrecord are drawn.  
        the mRNA and CDS are coupled by using a ``qualifier['note']`` containing 
        ``'transcript_id=XXXXXXXXXXXX'`` this work for ensembl-derived 
        seqrecord, but can be impoved.
        '''
        # get gene info
        for feat in self.seqrecord.features:
            if (feat.type=='gene'):
                if gene_code:
                    if gene_code==feat.qualifiers['gene'][0]:
                        Lgene_name=[gene_code]
                        if feat.qualifiers.has_key('locus'):
                            Lgene_name.extend(feat.qualifiers['locus'])
                        if feat.qualifiers.has_key('note'):
                            Lgene_name.extend(feat.qualifiers['note'])
                        gene_name=' - '.join(Lgene_name)
                        break
                else:
                    gene_code=feat.qualifiers['gene'][0]
                    Lgene_name=[gene_code]
                    if feat.qualifiers.has_key('locus_tag'):
                        Lgene_name.extend(feat.qualifiers['locus_tag'])
                    if feat.qualifiers.has_key('note'):
                        Lgene_name.extend(feat.qualifiers['note'])
                    gene_name=' - '.join(Lgene_name)
                    break
        if not gene_code:
            raise ValueError("Cannot draw gene structure, without a gene code")





        DFeature2Draw={'Gene':{'gene':[],'exons':[]},
                       'mRNA and CDS':[],
                       'Non coding RNA':[],
                       'Other':[],
                       'Feature not linked to gene':[],
                       }

        '''couple mRNA and CDS, and collect exons, misc_RNA and other features '''

        mRNA_codes=[]
        for feature in self.seqrecord.features:
            if feature.qualifiers.has_key('gene'):
                if  gene_code in feature.qualifiers['gene']:
                    if feature.type =='mRNA':
                        for note in feature.qualifiers['note']:#IMPROVE HERE
                            if 'transcript_id' in note:
                                transcript_id=note
                                mRNA_codes.append(note)
                        if feature.type =='mRNA':
                            for CDS in self.seqrecord.features:
                                if (CDS.type =='CDS') and (transcript_id in  CDS.qualifiers['note']):
                                    break
                            DFeature2Draw['mRNA and CDS'].append([feature,CDS])
                    elif feature.type !='CDS':
                        if feature.type =='gene':
                            DFeature2Draw['Gene']['gene'].append(feature)
                        elif feature.type =='exon':
                            DFeature2Draw['Gene']['exons'].append(feature)
                        elif feature.type =='misc_RNA':
                            DFeature2Draw['Non coding RNA'].append(feature)
                        else:
                            DFeature2Draw['Other'].append(feature)
            else:
                    DFeature2Draw['Feature not linked to gene'].append(feature)

        '''DRAW GENE'''
        grfeat=feature.GeneSeqFeature(DFeature2Draw['Gene']['gene'][0],DFeature2Draw['Gene']['exons'],type='Gene',name=gene_name,height=1,cm_value=0.8,alpha=0.5,cm='Purples',linewidth=1)
        self.panel.add_feature(grfeat,feature_class=1)
        '''DRAW CODING RNA'''
        c=0.
        for mRNA,CDS in DFeature2Draw['mRNA and CDS']:
            c+=1
            try:
                mRNA_name=', '.join(mRNA.qualifiers['note']).replace('_', ' ')+' - '+CDS.qualifiers['protein_id'][0]
            except:
                try:
                    mRNA_name=', '.join(mRNA.qualifiers['note']).replace('_', ' ')
                except:
                    mRNA_name='Transcript '+str(int(c))
            grfeat=features.CoupledmRNAandCDS(mRNA,CDS,alpha=0.75,type='mRNA',cm_value=c/len(DFeature2Draw['mRNA and CDS']),cm='autumn',name=mRNA_name)
            self.panel.add_feature(grfeat,feature_class=2)
        '''DRAW NON CODING RNA'''
        c=0.
        for RNA in DFeature2Draw['Non coding RNA']:
            c+=1
            try:
                RNA_name=', '.join(RNA.qualifiers['note']).replace('_', ' ')
            except:
                RNA_name='Non coding transcript '+str(int(c))
            grfeat=features.SegmentedSeqFeature(RNA,alpha=0.75,type='Non coding RNA',cm_value=c/len(DFeature2Draw['Non coding RNA']),cm='winter',name=RNA_name)
            self.panel.add_feature(grfeat,feature_class=3)
        '''DRAW OTHER FEATURES'''
        c=0.
        for feat in DFeature2Draw['Other']:
            c+=1
            try:
                feat_name=', '.join(mRNA.qualifiers['hit_name']).replace('_', ' ')
            except:
                feat_name=''
            grfeat=features.GenericSeqFeature(RNA,alpha=0.75,type='Other',cm_value=c/len(DFeature2Draw['Other']),name=feat_name)
            self.panel.add_feature(grfeat,feature_class=4)



def SliceSeqRec(seqrec,start=None,end=None,include_feature=[],exclude_feature=[],feature_rename_dict={}):
    ''' 
    Adapted from default biopython seqrecord slicing, features are retained 
    even if they are not completely included in the slice.

        ======================= ================================================
        Key                      Description
        ======================= ================================================
        include_feature         lists of feature to include from seqrecord, 
                                if empty all features are included
        
        exclude_feature         lists of feature to exclude from seqrecord, 
                                if empty no feature is excluded. 
                                `exclude_feature` overrides `include_feature`
        
        feature_rename_dict     key value dict used to change the `feature.type`
                                attributes to be correctly visualized in the 
                                drawer. 
                                eg.: 
                                ``feature_rename_dict={'Seg':'Low complexity'}`` 
        ======================= ================================================
    '''
    if seqrec.seq is None :
        raise ValueError("If the sequence is None, we cannot slice it.")
    parent_length = len(seqrec)
    '''Create new seqrecord'''
    answer = SeqRecord(seqrec.seq[start:end],
                       id=seqrec.id,
                       name=seqrec.name,
                       description=seqrec.description+' | Sliced:'+str(start)+'..'+str(end))
    '''Select relevant features, add them with shifted locations '''
    if start is None :
        start = 0
    if end is None :
        end = -1
    if (start < 0 or end < 0) and parent_length == 0 :
        raise ValueError("Cannot support negative indices without the sequence length")
    if start < 0 :
        start = parent_length - start
    if end < 0  :
        end  = parent_length - end + 1
    for f in seqrec.features :
        if (start <= f.location.end.position)  and (f.location.start.position < end):
            if (not include_feature) or f.type in include_feature:
                if f.type not in exclude_feature:
                    if f.type in feature_rename_dict:
                        f.type=feature_rename_dict[f.type]
                    answer.features.append(f._shift(-start))
    for key, value in seqrec.letter_annotations.iteritems() :
        answer._per_letter_annotations[key] = value[start:end]
    return answer

