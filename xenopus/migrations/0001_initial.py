# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('unigene_id', models.CharField(max_length=20)),
                ('protein_id', models.CharField(max_length=20)),
                ('gene_id', models.CharField(max_length=20)),
                ('protein_length', models.PositiveIntegerField(default=0)),
                ('gene_length', models.PositiveIntegerField(default=0)),
                ('protein_name', models.CharField(max_length=200)),
                ('gene_name', models.CharField(max_length=200)),
                ('organism', models.CharField(max_length=200)),
                ('gene_type', models.CharField(max_length=20)),
                ('protein_similarity', models.FloatField(default=0, max_length=5)),
                ('protein_classification', models.CharField(max_length=10)),
                ('protein_tm_helices', models.PositiveIntegerField(default=0)),
                ('protein_predicted_segments', models.TextField()),
                ('gene_start_nucleotide', models.PositiveIntegerField(default=0)),
                ('gene_end_nucleotide', models.PositiveIntegerField(default=0)),
                ('gene_start_codon', models.CharField(max_length=3)),
                ('gene_stop_codon', models.CharField(max_length=3)),
                ('protein_start_aa', models.CharField(max_length=1)),
                ('protein_stop_aa', models.CharField(max_length=1)),
                ('gene_sequence', models.TextField()),
                ('protein_sequence', models.TextField()),
            ],
        ),
    ]
