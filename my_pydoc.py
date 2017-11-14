import django
import pydoc
import os

os.environ['DJANGO_SETTINGS_MODULE'] = 'iac.settings'
django.setup()
pydoc.cli()
