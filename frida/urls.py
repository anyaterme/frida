from django.conf.urls import include, url
from frida import views

urlpatterns = [
	url(r'^$', views.index, name='index'),
	url(r'^calculate/$', views.calculate, name='calculate'),
	url(r'^doc/params1/$', views.docs_p1, name='docs-p1'),
]
