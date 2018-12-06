# -*- encoding: utf-8 -*-
from django import template
from django.conf import settings
from django.utils.safestring import mark_safe
from django.urls import reverse
import datetime
import numpy as np

register = template.Library()

@register.filter
def percentage(value, decimals=0):
	try:
		value = np.round(value*100, decimals)
		return value
	except:
		return -1

@register.filter
def astro_round(quantity, decimals=0):
	try:
		myfmt = "{0.value:0.0%df} {0.unit}" % decimals
		return myfmt.format(quantity)
	except:
		myfmt = "{0.value:0.0%df} {0.unit}" % decimals
		return myfmt.format(quantity[0])

@register.filter
def astro_science(quantity, decimals=0):
	try:
		myfmt = "{0.value:0.0%de} {0.unit}" % decimals
		return myfmt.format(quantity)
	except:
		try:
			myfmt = "{0.value:0.0%E} {0.unit}" 
			return myfmt.format(quantity[0])
		except:
			return "0"

@register.simple_tag
def astro_scien_to(quantity, decimals=4, unit=None):
	try:
		quantity = quantity.to(unit)
		myfmt = "{0.value:0.0%de} {0.unit}" % decimals
		return myfmt.format(quantity)
	except Exception as e:
		print (e)
		try:
			myfmt = "{0.value:0.0%de} {0.unit}" % decimals
			return myfmt.format(quantity[0])
		except Exception as e:
			print (e)
			return "0"
		myfmt = "{0.value:0.0%de} {0.unit}" % decimals
		return myfmt.format(quantity[0])

@register.simple_tag
def astro_round_to(quantity, decimals=0, unit=None):
	if unit is not None:
		quantity = quantity.to(unit)
	try:
		myfmt = "{0.value:0.0%df} {0.unit}" % decimals
		return myfmt.format(quantity)
	except:
		myfmt = "{0.value:0.0%df} {0.unit}" % decimals
		return myfmt.format(quantity[0])

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



