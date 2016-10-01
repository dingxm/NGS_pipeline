#!/usr/bin/env python
"""
"""
import mechanize
import cookielib
#import html2text
import re
import time
import sys

hgmd_login_url = "http://www.hgmd.cf.ac.uk/docs/login.html"
email_address = "dingxm@ucla.edu"
password = "HGMD189066"
if email_address == "":
    sys.exit("define your email address and password")



def initialize_browser():
    """
    Initialize a Browser object

    thanks to http://stockrt.github.com/p/emulating-a-browser-in-python-with-mechanize/
    """

    br = mechanize.Browser()
    # Cookie Jar
    cj = cookielib.LWPCookieJar()
    br.set_cookiejar(cj)

    # Browser options
    br.set_handle_equiv(True)
#    br.set_handle_gzip(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)

    # Follows refresh 0 but not hangs on refresh > 0
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)

    # Want debugging messages?
#    br.set_debug_http(True)
#    br.set_debug_redirects(True)
#    br.set_debug_responses(True)

    # User-Agent (this is cheating, ok?)
    br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
    br.addheaders.append(('email', email_address))

    return br


def login_hgmd(br):
    """
    login to HGMD

    After calling this function, you will be able to search in HGMD programmatically
    """
    response = br.open(hgmd_login_url)
    html = response.read()

    # print response to STDOUT for debugging purposes
    # the html2text library is used for formatting the output in a more readable form
    print html

    # print all the forms in the current page
    print [f for f in br.forms()]

    # select login form
    br.select_form(nr=0)
    print br.form

    # print all controls in the current form, for debugging purposes
    print [c.name for c in br.form.controls]

    # set username and password
    br.form['email'] = email_address
    br.form['password'] = password

    # submit form
    response_form = br.submit()

    # Now, you should have successfully logged in. The contents of the page will be changed. Check the contents of br.read()
    html_response = response_form.read()
    print html_response

    # Then, you should complete this on your own. I suggest you to br.open("http://www.hgmd.cf.ac.uk/ac/index.php"), select the Search form, and submit a query again

    # wait 2 seconds to not overload the server
    time.sleep(2)

    return br


#def pretty_print_page(br):
#    print html2text.html2text(br.response().read())



if __name__ == '__main__':
    br = initialize_browser()
#    resp = browse_dbcline(br, genes = ['GCS1'])
    br = login_hgmd(br)
