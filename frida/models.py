from django.db import models
from django.utils.translation import ugettext_lazy as _

class Test(models.Model):
	name = models.CharField(max_length=200, verbose_name=_("Name"))

	def __unicode__(self):
		return self.name

	class Meta:
		verbose_name = _('Test')
		verbose_name_plural = _('Tests')
