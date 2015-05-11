# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('xenopus', '0002_auto_20150511_0107'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='protein_similarity',
            field=models.FloatField(default=0, max_length=5),
        ),
    ]
