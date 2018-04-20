# -*- encoding: utf-8 -*-
from django import template
from django.conf import settings
from django.utils.safestring import mark_safe
from django.core.urlresolvers import reverse
import datetime

register = template.Library()

@register.filter
def percentage(value, decimals=0):
	try:
		value = round(value*100, decimals)
		return value
	except:
		return "ERROR percentage"

@register.filter
def list_value(dict, index):
	return dict[index]

@register.filter
def ellipsis_text(text, size):
	if (len(text) > size):
		return("%s..." % text[:size])
	else:
		return text

@register.filter
def astro_unit_value(item, unit=None):
	if unit is None:
		unit =item.unit
	return item.to(unit).value

@register.filter
def astro_unit(item, unit):
	return item.to(unit)


@register.simple_tag
def ajax_send_form_reverse(action, form, div):
	url = reverse(action)
	return ("$.ajax({url:'%s',type:'POST',data:$('#%s').serialize(),success:function(data){$('#%s').html(data);gotoelement('%s');}});" % (url, form, div, div))



