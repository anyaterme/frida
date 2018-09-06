from django.conf.urls import include, url
from frida import views

urlpatterns = [
	url(r'^$', views.index, name='index'),
	url(r'^calculate/$', views.calculate, name='calculate'),
	url(r'^doc/params1/$', views.docs_p1, name='docs-p1'),
	url(r'^doc/params2/$', views.docs_p2, name='docs-p2'),
	url(r'^doc/params3/$', views.docs_p3, name='docs-p3'),
]
