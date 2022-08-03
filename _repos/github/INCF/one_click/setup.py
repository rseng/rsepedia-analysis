from setuptools import setup
 
APP = ['qtpush.py']
OPTIONS = {'argv_emulation': True, 
           'includes': ['sip', 'PyQt4'], 
           'excludes': ['PyQt4.uic.port_v3'], 
           'iconfile': 'mac_icon.icns'}
 
setup(
    app=APP,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
    data_files = [ ('', ['agreement_dialog.ui', 
                         'exit_notice_dialog.ui', 
                         'file_loader_dialog.ui', 
                         'info_dialog.ui', 
                         'login_dialog.ui', 
                         'quit_dialog.ui', 
                         'upload_dialog.ui']) ],
)
