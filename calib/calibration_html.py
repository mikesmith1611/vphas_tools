import numpy as np
from astropy.table import Table
import os

reference = 'panstarrs'

def create_template_html():
    html_str = """
<html>

<head>

    <script type="text/javascript" src="sorttable.js"></script>
    <LINK href="style.css" rel="stylesheet" type="text/css">

</head>

<body>

<h1> VPHAS+ calibration shifts</h1>
<a href="tools/calib/tables/apass/shifts_apass_all.fits">Full Table (fits)</a> &nbsp;&nbsp;&nbsp;&nbsp;
<a href="tools/calib/plots/">All plots</a> &nbsp;&nbsp;&nbsp;&nbsp;

<br> <br>

<table class="sortable" border="1">
<tr>
<th>Field</th>
<th>Concat</th>
<th>l</th>
<th>b</th>
<th>u shift</th>
<th>g shift</th>
<th>r1 shift</th>
<th>r2 shift</th>
<th>i shift</th>
<th>u concat shift </th>
<th>g concat shift</th>
<th>r1 concat shift</th>
<th>r2 concat shift</th>
<th>i concat shift</th>
<th>Pointing plots</th>
<th>Concat plots</th>

</tr>
        """
    html_file = open('/home/msmith/public_html/vphasshifts_{0}.html'.format(reference), 'w')
    html_file.write(html_str)
    html_file.close()

def add_row(field, data, i, makethumb=False, make_html_image=True):

    mask = (data['field'] == f)
    data2 = data
    data = data[mask]

    paramimage = 'tools/calib/plots/{0}/vphas_{1}_{0}.png'.format(reference, field)
    concat_img = 'tools/calib/plots/{0}/vphas_{1}_concat_calib.png'.format(reference, field)
    paramimagesm = paramimage[:-4] + '_s.png'

    if make_html_image:
        paramimage = make_image_html(i, paramimage, data2)
        concat_img = make_image_html2(i, concat_img, data2)
        print paramimage
    if makethumb:
        os.system('convert ' + paramimage + ' -resize 100x175 ' + paramimagesm)
        os.system('convert ' + paramimage[:-4] + '_s.pdf ' + paramimage[:-4] +
                  '_s.png')

    l = "%.5f" % data['GAL_LONG'][0]
    b = "%.5f" % data['GAL_LAT'][0]

    u = "%.2f" % data['u'][0]
    g = "%.2f" % data['g'][0]
    r = "%.2f" % data['r1'][0]
    r2 = "%.2f" % data['r2'][0]
    i = "%.2f" % data['i'][0]
    u_t = "%.2f" % data['u_c'][0]
    g_t = "%.2f" % data['g_c'][0]
    r_t = "%.2f" % data['r1_c'][0]
    r2_t = "%.2f" % data['r2_c'][0]
    i_t = "%.2f" % data['i_c'][0]
    concat = str(data['Concat'][0])

    html_str = r'<tr>' + '\n'
    html_str += r'<td>' + str(field) + '</td>' + '\n'
    html_str += r'<td>' + concat + '</td>' + '\n'

    #html_str += r'<td>' + str(ugrdate) + '</td>' + '\n'
    #html_str += r'<td>' + str(haridate) + '</td>' + '\n'

    html_str += r'<td>' + l + '</td>' + '\n'
    html_str += r'<td>' + b + '</td>' + '\n'

    html_str += r'<td>' + u + '</td>' + '\n'
    html_str += r'<td>' + g + '</td>' + '\n'
    html_str += r'<td>' + r + '</td>' + '\n'
    html_str += r'<td>' + r2 + '</td>' + '\n'
    html_str += r'<td>' + i + '</td>' + '\n'
    html_str += r'<td>' + u_t + '</td>' + '\n'
    html_str += r'<td>' + g_t + '</td>' + '\n'
    html_str += r'<td>' + r_t + '</td>' + '\n'
    html_str += r'<td>' + r2_t + '</td>' + '\n'
    html_str += r'<td>' + i_t + '</td>' + '\n'


    html_str += r'<td><a href="' + paramimage + '" target="_blank">Pointing plots</a>' + '\n'  #  '<img src="' + paramimagesm + '"width="100"></a></td>' + '\n'
    html_str += r'<td><a href="' + concat_img + '" target="_blank">Concat plots</a>' + '\n'
    html_str += r'</tr>' + '\n'

    html_file = open('/home/msmith/public_html/vphasshifts_{0}.html'.format(reference), 'a')
    html_file.write(html_str)
    html_file.close()


def end_html():

    html_str = """
                </table>
                </body>
                </html> \n
            """
    html_file = open('/home/msmith/public_html/vphasshifts_{0}.html'.format(reference), 'a')
    html_file.write(html_str)
    html_file.close()


def make_image_html(i, image, data):
    html_str = "<html> \n <head> \n"

    html_str += r"""

<script language="javascript">
function sayHello()
{
alert("Flagged !");
return true;
}
</script>
<LINK href="btn.css" rel="stylesheet" type="text/css">
</head>
<body>
"""
    image2 = image
    html_str += '<embed src="../../' + image2 + '" style="float: left; width: 60%;" />'
    html_str += '<embed src="../../tools/calib/plots/{1}/spatial_map_{0}.png"'.format(data['field'][i], reference)
    html_str += ' style="float: left; width: 40%;" /> \n'
    html_str += '<p style="clear: both;">'
    html_str += '<br><p align="center"> '
    if i == 0:
        html_str += ('<a href="vphas_' + data['field'][-1] + '_' + reference +
                     '.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][i + 1] + '_' + reference +
                     '.html' + '">Next</a>')

    elif i != len(data) - 1:

        html_str += ('<a href="vphas_' + data['field'][i - 1] + '_' + reference +
                     '.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][i + 1] + '_' + reference +
                     '.html' + '">Next</a>')

    else:

        html_str += ('<a href="vphas_' + data['field'][i - 1] + '_' + reference +
                     '.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][0] + '_' + reference +
                     '.html' + '">Next</a>')

    html_str += '</p>'
    html_str += '<a href="vphas_' + data['field'][i] + '_concat_calib.html" > Concat Plots</a>'

    html_str += '<form action="save.php?field=' + data['field'][i] +'" method="post" onSubmit="return sayHello()">'
    html_str += '<input  class="btn" type="submit" Value="Flag !" align="right">'
    html_str += '</form>'

    html_str += "</body> \n </html>"
    img = image.split('/')[-1].split('.')[0]
    html_fname = '/home/msmith/public_html/calibplots/{0}/'.format(reference) + img + '.html'
    html_file = open(html_fname, 'w')
    html_file.write(html_str)
    html_file.close()
    return 'calibplots/{0}/'.format(reference) + img + '.html'

def make_image_html2(i, image, data):
    html_str = "<html> \n <head> \n"

    html_str += r"""

<script language="javascript">
function sayHello()
{
alert("Flagged !");
return true;
}
</script>
<LINK href="btn.css" rel="stylesheet" type="text/css">
</head>
<body>
"""
    image2 = image
    html_str += '<div style="float:left;"><embed src="../../' + image2 + '" height="90%"/> \n'
    html_str +=  '<div style="float:right;"><table style="width:50%" border="1"> \n <tr>'
    html_str += '<th>u</td> \n'
    html_str += '<th>g</td> \n'
    html_str += '<th>r1</td> \n'
    html_str += '<th>r2</td> \n'
    html_str += '<th>i</td> \n </tr>\n <tr>'

    html_str += '<td>{0:.2f}</td> \n'.format(data['u_c'][i])
    html_str += '<td>{0:.2f}</td> \n'.format(data['g_c'][i])
    html_str += '<td>{0:.2f}</td> \n'.format(data['r1_c'][i])
    html_str += '<td>{0:.2f}</td> \n'.format(data['r2_c'][i])
    html_str += '<td>{0:.2f}</td> \n </tr></table></div>'.format(data['i_c'][i])
    #html_str += '</p>'
    #html_str += '<p>'
    html_str += '<br><p align="center"> '
    if i == 0:
        html_str += ('<a href="vphas_' + data['field'][-1] +
                     '_concat_calib.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][i + 1] +
                     '_concat_calib.html' + '">Next</a>')

    elif i != len(data) - 1:

        html_str += ('<a href="vphas_' + data['field'][i - 1] +
                     '_concat_calib.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][i + 1] +
                     '_concat_calib.html' + '">Next</a> \n')

    else:

        html_str += ('<a href="vphas_' + data['field'][i - 1] +
                     '_concat_calib.html' + '">Previous | </a> ')

        html_str += ('<a href="vphas_' + data['field'][0] +
                     '_concat_calib.html' + '">Next</a> \n')
    html_str += '</p>'
    html_str += '<a href="vphas_' + data['field'][i] + '.html" > Single Calib Plots</a>'
    html_str += '<form action="save.php?field=' + data['field'][i] +'" method="post" onSubmit="return sayHello()">'
    html_str += '<input  class="btn" type="submit" Value="Flag !" align="right">'
    html_str += '</form>'
    html_str += "</body> \n </html>"
    img = image.split('/')[-1].split('.')[0]
    html_fname = '/home/msmith/public_html/calibplots/{0}/'.format(reference) + img + '.html'
    html_file = open(html_fname, 'w')
    html_file.write(html_str)
    html_file.close()
    return 'calibplots/{0}/'.format(reference) + img + '.html'


create_template_html()
"""
for i, j in enumerate(data['ID']):
    os.system('mv ' + '../data/OBstars/good/' + j.strip() + '*.png ' +
              '../data/OBstars/emission/')
"""
data = Table.read('tables/{0}/shifts_{0}_all.fits'.format(reference))
data.sort('field')
for i, f in enumerate(data['field']):
    add_row(f, data, i, makethumb=False, make_html_image=True)

end_html()

"""
files = glob.glob('../data/shifts/vphas_apass_????.pdf')
it = 0
pb = vphas.progress_bar(len(files))
for i, j in enumerate(files):
    os.system('convert -density 150 ' + j + ' -quality 100 ' + j.strip('pdf') + 'png')
    it += 1
    pb.update(it)
"""
