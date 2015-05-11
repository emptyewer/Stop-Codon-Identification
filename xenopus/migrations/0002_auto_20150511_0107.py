# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('xenopus', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='gene',
            name='gene_organism',
        ),
        migrations.AlterField(
            model_name='gene',
            name='gene_end_nucleotide',
            field=models.PositiveIntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='gene',
            name='gene_length',
            field=models.PositiveIntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='gene',
            name='gene_start_nucleotide',
            field=models.PositiveIntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='gene',
            name='protein_length',
            field=models.PositiveIntegerField(default=0),
        ),
    ]
