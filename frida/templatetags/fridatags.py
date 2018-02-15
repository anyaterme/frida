# -*- encoding: utf-8 -*-
from django import template
from django.conf import settings
from django.utils.safestring import mark_safe
import datetime

register = template.Library()

@register.filter
def list_value(dict, index):
	return dict[index]

@register.filter
def astro_unit_value(item, unit):
	return item.to(unit).value
