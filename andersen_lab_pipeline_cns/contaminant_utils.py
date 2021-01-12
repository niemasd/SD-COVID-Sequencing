import re  # regex
import numpy as np  # linear algebra
import pandas as pd  # data processing
from path import Path  # file I/O
import plotly  # viz
import plotly.figure_factory as ff  # viz
import plotly.express as px  # viz
import plotly.graph_objects as go  #viz
from plotly.subplots import make_subplots
from jinja2 import Environment, FileSystemLoader  # html template engine




def generate_html(general_hmap, cont_hmap, cont_table, num_samples, num_conts, expt_name):
    # express plots in html and JS
    general_hmap = plotly.offline.plot(general_hmap, include_plotlyjs=False, output_type='div')
    cont_hmap = plotly.offline.plot(cont_hmap, include_plotlyjs=False, output_type='div')
    # express contaminants table in html
    cont_table = cont_table.to_html()
    # dir containing our template
    file_loader = FileSystemLoader('templates')
    # load the environment
    env = Environment(loader=file_loader)
    # load the template
    template = env.get_template('contamination.html')
    # render data in our template format
    html_output = template.render(general_hmap=general_hmap, cont_hmap=cont_hmap, 
                                  cont_table=cont_table, num_samples=num_samples,
                                  num_conts=num_conts, expt_name=expt_name)
    return html_output


def save_html(html_output: str, filename: str):
    with open(filename, 'w') as f:
        f.write(html_output)
        
        
def generate_heatmap(data, x, y, normalize=True):
    counts = data[:, :-1]
    if normalize:
        # normalize counts per sample (to address visualization issue) - IGNORE IF USING LOG
        summed_counts = counts.sum(axis=1)[:, np.newaxis]
        # normalize counts [OPTIONAL]
        counts_normed = counts / np.where(summed_counts > 0, summed_counts, 1)
    # grab column with contaminant flags
    flags = data[:, -1][:, np.newaxis]
    # generate plot
    fig = make_subplots(2,1, shared_xaxes=True)

    fig.add_trace(
     go.Heatmap(z = counts_normed.T, x = x, y = y, hovertext=counts.T, hovertemplate =
                '<b>Sample</b>: %{y}<br>'+
                '<b>Barcode</b>: %{x}<br>'+
                '<b>Count</b>: %{hovertext}<br>'+
                '<b>Normed Count</b>: %{z:.2f}',
                hoverinfo='all', coloraxis = "coloraxis"), 2,1)

    fig.add_trace(go.Heatmap(z = flags.T, x=x, y=['Contamination'], coloraxis = "coloraxis"),1,1)
    fig.update_layout(coloraxis = {'colorscale':'viridis'}, height=800,
        yaxis=dict(
            domain=[0.97, 1]
        ),
        yaxis2=dict(
            domain=[0, 0.95]
        )
    )
    return fig


def get_heatmap_data(ans: pd.DataFrame):
    """Sorts the values for the heatmap to make it more intuitive and readable
    Expects a dataframe `ans` containing processed (clean) barcode read counts data
    Expects a dataframe `flag`"""
    # generate heatmap matrix of (logged) read counts per sample per paired read
    hmap = (ans.pivot_table(index=["sample"], columns=["paired_read"], values="paired_read_count")
               .replace([np.inf, -np.inf], np.nan)
               .fillna(0))
    # drop unknown-unknown reads only
    hmap = hmap.drop(columns='unknown-unknown')
    # sort values to make plot more intuitive
    hmap['max_idx'] = hmap.apply(lambda x: hmap.columns.tolist().index(x.idxmax()), axis=1)
    hmap = hmap.sort_values('max_idx').drop(columns=['max_idx'])
    # grab read counts
    counts = hmap.values
    # prepare data for identifying potential contaminants
    flag = (ans.groupby('sample')
               .agg(uniq_forward_bcodes = ('forward_barcode', get_unique_barcodes),
                    uniq_reverse_bcodes = ('reverse_barcode', get_unique_barcodes)))

    # create boolean column that identifies potential contamination
    flag['contamination'] = flag.apply(is_contaminant, axis=1)
    # merge with original data to include the contamination flags
    contaminants_flag = (hmap.join(flag, how='inner')['contamination']
                             .apply(lambda x: np.where(x==True, 1.0, 0.0)))
    # add contaminant flag column to the read counts
    data = np.hstack((counts, contaminants_flag[:, np.newaxis]))
    # list of all sample IDs
    x = hmap.index.values
    # list of all paired reads and an extra flag column for contamination
    y = hmap.columns.tolist()
    return hmap, data, x, y


def get_contaminated_data(data: np.array, hmap: pd.DataFrame):
    # choosing only samples with potential contaminants
    cont_data = data[data[:, -1] !=0][:, :-1]
    # list of samples with potential contamination only
    cont_y = hmap.columns.values[~np.all(cont_data == 0, axis=0)].tolist()
    # data (read counts) for contaminated samples only
    cont_data = cont_data[:, ~np.all(cont_data == 0, axis=0)]
    # list of paired reads with potential contamination only
    cont_x = hmap.index.values[data[:, -1] !=0].tolist()
    return np.hstack((cont_data, np.ones((cont_data.shape[0], 1)))), cont_x, cont_y
        
# def generate_heatmap(data, x, y):
#     heatmap = go.Heatmap(z=data.T, x=x, y=y)
#     plot = [heatmap]
#     fig = go.Figure(data = plot)
#     return plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')



def generate_table(ans, contaminated_samples):
    return (ans.loc[(ans['sample'].isin(contaminated_samples)) & (ans['paired_read_count'] > 0)]
                 .set_index(['sample', 'paired_read'])
                 .sort_index()[['paired_read_count']])


def load_all_data(input_pths: list):
    """Load barcode read data from the given list of sample file paths"""
    ans = pd.DataFrame()
    for sample_pth in input_pths:
        try:
            df = prepare_data(*load_data(sample_pth))
        except: continue
        df = (df.reset_index()#[['forward_barcode', 'reverse_barcode']]
                .drop_duplicates())
        df['sample'] = sample_pth.basename().split('_')[0]
        ans = pd.concat([ans, df], axis=0)
    return ans


def load_data(input_file: Path):
    with open(input_file, 'r') as f:
        input_data = f.readlines()
    header = re.split('\t|\n', input_data[0])[:-1]
    data = []
    for row in input_data[1:]:
        data.append(re.split('\t|\n', row)[:-1])
    df = pd.DataFrame(data=data, columns=header)
    return df, header


def prepare_data(df: pd.DataFrame, header: list):
    # convert read counts to integer
    df[header[-1]] = df[header[-1]].astype(int)
    # consolidate barcodes and their reverse complements 
    df = df.apply(_consolidate_reverse_complements, axis=1)
    # merge barcodes and their rcs and sum their read counts
    df = (df.groupby(['forward_barcode', 'reverse_barcode'])
            .agg({'paired_read_count': 'sum'})
            .reset_index())
    # create df with all possible paired combos of forward and reverse barcodes 
    forward_bcodes = df[['forward_barcode']].drop_duplicates()
    reverse_bcodes = df[['reverse_barcode']].drop_duplicates()
    forward_bcodes['key'] = 0
    reverse_bcodes['key'] = 0
    all_pairs = forward_bcodes.merge(reverse_bcodes, how='outer', on='key').drop(columns='key')
    all_pairs = (all_pairs.merge(df, on=['forward_barcode', 'reverse_barcode'], how='left')
                          .fillna(0)
                          .set_index(['forward_barcode', 'reverse_barcode']))
    return all_pairs

def _consolidate_reverse_complements(x: str):
    x['forward_barcode'] = x['forward_barcode'].split('_')[0]
    x['reverse_barcode'] = x['reverse_barcode'].split('_')[0]
    return x


def get_unique_barcodes(x):
    x = set(x.unique())
    x.discard('unknown')
    return x


def is_contaminant(x):
    if len(x['uniq_forward_bcodes']) > 1 or len(x['uniq_reverse_bcodes']) > 1:
        return True
    return False