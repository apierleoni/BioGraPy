import unittest
import os
import tempfile
import Image
from Bio import SeqIO
from plone4bio.graphics import Panel
from plone4bio.graphics import features

# inspired by BioPerl Tutorial from http://bioperl.org/wiki/HOWTO:Graphics

class TestTutorial(unittest.TestCase):
    def setUp(self):
        # hit           score   start   end
        self.data1 = [
            ('hsHOX3', '381', 2, 200),
            ('scHOX3', '210', 2, 210),
            ('xlHOX3', '800', 2, 200),
            ('hsHOX2', '1000', 380, 921),
            ('scHOX2', '812', 402, 972),
            ('xlHOX2', '1200', 400, 970),
            ('BUM', '400', 300, 620),
            ('PRES1', '127', 310, 700),
        ]

    def test_example1(self):
        """
        # This is code example 1 in the Graphics-HOWTO
        use strict;
        use Bio::Graphics;
        use Bio::SeqFeature::Generic;
        my $panel = Bio::Graphics::Panel->new(-length => 1000, -width  => 800);
        my $track = $panel->add_track(-glyph => 'generic', -label => 1);
        while (<>) { # read blast file
            chomp;
            next if /^\#/;  # ignore comments
            my($name,$score,$start,$end) = split /\t+/;
            my $feature = Bio::SeqFeature::Generic->new(
                                      -display_name => $name,
                                      -score        => $score,
                                      -start        => $start,
                                      -end          => $end
                                   );
            $track->add_feature($feature);
        }
        print $panel->png;
        """
        # TODO: drawer vs. panel+track
        panel = Panel(xmin = 1, xmax = 1000)
        for data in self.data1:
            feature = features.Simple(name=data[0], score=data[1], start=data[2], end=data[3])
            panel.add_feature(feature)
        fh = tempfile.TemporaryFile()
        panel.save(fh, format='png')
        fh.seek(0)
        img = Image.open(fh)
        #self.assertEqual(img.size, (1000, 300))
        self.assertEqual(img.format, 'PNG')
        # panel.save('/tmp/xxx.png')
        # os.system('firefox /tmp/xxx.png')

    # TODO: implements pad_left/right
    def test_example2(self):
        """
        Example 2. Rendering the blast hit file with scores and scale

        #!/usr/bin/perl
        # This is code example 2 in the Graphics-HOWTO
        use strict;
        use lib '/home/lstein/projects/bioperl-live';
        use Bio::Graphics;
        use Bio::SeqFeature::Generic;
        my $panel = Bio::Graphics::Panel->new(
                               -length    => 1000,
                               -width     => 800,
                               -pad_left  => 10,
                               -pad_right => 10,
                               );
        my $full_length = Bio::SeqFeature::Generic->new(
                               -start => 1,
                               -end   => 1000,
                               );
        $panel->add_track($full_length,
                               -glyph   => 'arrow',
                               -tick    => 2,
                               -fgcolor => 'black',
                               -double  => 1,
                               );
        my $track = $panel->add_track(
                               -glyph     => 'graded_segments',
                               -label     => 1,
                               -bgcolor   => 'blue',
                               -min_score => 0,
                               -max_score => 1000,
                               );
        while (<>) { # read blast file
            chomp;
            next if /^\#/;  # ignore comments
            my($name,$score,$start,$end) = split /\t+/;
            my $feature = Bio::SeqFeature::Generic->new(
                               -display_name => $name,
                               -score        => $score,
                               -start        => $start,
                               -end          => $end
                               );
            $track->add_feature($feature);
        }
        print $panel->png;
        """
        # TODO
        pass

def test_suite():
    return unittest.makeSuite(TestTutorial)

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')

