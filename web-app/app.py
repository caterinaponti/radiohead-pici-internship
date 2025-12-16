from flask import Flask, render_template, request
import pandas as pd
import altair as alt
alt.data_transformers.disable_max_rows()
#alt.data_transformers.enable("vegafusion")
import matplotlib as mpl
import os

mpl.use("agg")

app = Flask(__name__)

# Paths
REDCAP_PATH = 'datasets/'
BIOREPO_PATH = 'datasets/'

# Load data
df_clean = pd.read_csv(os.path.join(REDCAP_PATH, 'df_clean_ready_for_web_app.csv'))
biorepo_2024_df = pd.read_csv(os.path.join(BIOREPO_PATH, 'biorepo_2024_df_clean_REDCap_records_merged.csv'))
biorepo_2025_df = pd.read_csv(os.path.join(BIOREPO_PATH, 'biorepo_2025_df_clean_REDCap_records_merged.csv'))
biorepo_2024_df = biorepo_2024_df.drop(columns=['Unnamed: 0'])
iorepo_2025_df = biorepo_2025_df.drop(columns=['Unnamed: 0'])

pretreat_columns = [c for c in df_clean.columns if 'pretreat' in c]
days_disorders_columns = [c for c in df_clean.columns if 'days' in c]
irAE_columns = [c for c in df_clean.columns if 'irAE' in c or 'irae' in c]
for col in ['min_irae_days_first', 'Steroid_irAE_NA99', 'irae_n']:
    if col in irAE_columns:
        irAE_columns.remove(col)

general_medical_history_cols = [
    'subject_metastatic_disease', 'ever_smoked', 'prior_targeted_therapy',
    'prior_surgery', 'prior_radiation', 'prior_chemotherapy', 'prior_treated',
    'prior_car_t', 'respiratory disease', 'blood condition', 'cardiac/heart disease',
    'abnormal blood pressure', 'allergies', 'gastrointestinal problems',
    'musculoskeletal condition/problem', 'neurological problem',
    'endocrine condition/problem', 'diabetes', 'vision condition/problem',
    'sinus condition/problem', 'sleep', 'previous cancer',
    'dematological conditions', 'kidney/genitourinary disease', 'sleep apnea',
    'hearing condition/problem', 'infectious disease', 'raynauds',
    'antibiotics_yn', 'infection_hist_yn', 'pneumo_vax_hist', 
    'zoster_vax_hist', 'baseline_influenza_vax'
]

color_by_cols = [
    'stage', 'line_1_group', 'subject_dead', 'subject_dead_hospice',
    'subject_os_detail', 'subject_dead_hospice_censor',
    'clinical_observation_os_event', 'subject_os_detail_censor',
    'clinical_observation_pfs_event', 'probable_pfs_event', 'offstudy_reason',
    'infection_hist_post_abx', 'post_15_toxevent', 'subject_ethnicity',
    'line_1_therapy', 'subject_race', 'smoking_status', 'subject_sex',
    'min_irae_days_first', 'Steroid_irAE_NA99'
]



@app.route('/', methods=['GET', 'POST'])
def index(): 
    #select samples 
    selected_biorepository = request.form.get('biorepository', '2024')
    biorepo_df = get_biorepo_df(selected_biorepository)
    
    timepoints = sorted(biorepo_df["Barcode Timepoint"].dropna().unique())
    cancer_types = sorted(df_clean["meddra_disease_preferred_name"].dropna().unique())

    # Defaults
    selected_medical_history = request.form.get('medical_history', 'All')
    selected_pre_treat = request.form.get('pre_treat', 'All')
    selected_irAE = request.form.get('irae', 'All')
    x_axis_choice = request.form.get('x_axis', 'Cancer Type')
    color_by = request.form.get('color_by', 'None')
    if color_by not in color_by_cols:
        color_by = 'None'

    # Filter
    filtered_df = biorepo_df.copy()
    if selected_medical_history != "All":
        ids = df_clean[df_clean[selected_medical_history] == 1]["subject_id"].unique()
        filtered_df = filtered_df[filtered_df["Participant ID"].isin(ids)]
    if selected_pre_treat != "All":
        ids = df_clean[df_clean[selected_pre_treat] == 1]["subject_id"].unique()
        filtered_df = filtered_df[filtered_df["Participant ID"].isin(ids)]
    if selected_irAE != "All":
        ids = df_clean[df_clean[selected_irAE] == 1]["subject_id"].unique()
        filtered_df = filtered_df[filtered_df["Participant ID"].isin(ids)]


    # Merge
    merge_cols = ["subject_id", "meddra_disease_preferred_name"]
    if color_by != 'None':
        merge_cols.append(color_by)
    barplot_df = filtered_df.merge(
        df_clean[merge_cols],
        how="left",
        left_on="Participant ID",
        right_on="subject_id"
    )

    x_var = "meddra_disease_preferred_name" if x_axis_choice == "Cancer Type" else "Barcode Timepoint"
    barplot_df = barplot_df.dropna(subset=[x_var])

    # Bar Chart
    bar = alt.Chart(barplot_df).mark_bar().encode(
        x=alt.X(f'{x_var}:N', sort='-y'),
        y='count():Q',
        tooltip=[f"{x_var}:N", 'count():Q']
    )
    if color_by != "None" and color_by in barplot_df.columns:
        bar = bar.encode(color=alt.Color(f'{color_by}:N'))

    bar = bar.properties(width=600, height=400).interactive()
    bar_html = bar.to_html()

    return render_template(
        "index.html",
        selected_biorepository=selected_biorepository,
        bar_chart=bar_html,
        pretreat_columns=pretreat_columns,
        irAE_columns=irAE_columns,
        general_medical_history_cols=general_medical_history_cols,
        color_by_cols=color_by_cols + pretreat_columns + irAE_columns + general_medical_history_cols,
        x_axis_choice=x_axis_choice,
        selected_medical_history=selected_medical_history,
        selected_pre_treat=selected_pre_treat,
        selected_irAE=selected_irAE,
        color_by=color_by,
        filtered_df=filtered_df.drop(columns=['Unnamed: 0'], errors='ignore')
    )


def get_biorepo_df(selected_biorepository):
    if selected_biorepository == '2025':
        return biorepo_2025_df
    return biorepo_2024_df


@app.route('/material_type', methods=['POST'])
def submit_action():
    if request.method == 'POST':
        # Process the button click and any submitted data
        # For example:
        # data = request.form['input_field_name']

        selected_biorepository = request.form.get('biorepository', '2024')
        biorepo_df = get_biorepo_df(selected_biorepository)

        timepoints = sorted(biorepo_df["Barcode Timepoint"].dropna().unique())
        cancer_types = sorted(df_clean["meddra_disease_preferred_name"].dropna().unique())

        # Donut Plot with Altair 
        timepoint = request.form.get('timepoint', 'All')
        cancer = request.form.get('cancer', 'All')

        filtered_df2 = biorepo_df.copy()
        # Use Flask form data to determine filter2 value
        filter2 = request.form.get('patient_tubes', 'Patient Count')

        if timepoint != "All":
            filtered_df2 = filtered_df2[filtered_df2["Barcode Timepoint"] == timepoint]
        
        if cancer != "All":
            ids = df_clean[df_clean["meddra_disease_preferred_name"] == cancer]["subject_id"].unique()
            filtered_df2 = filtered_df2[filtered_df2["Participant ID"].isin(ids)]

        
        filtered_df2["Group"] = filtered_df2["Participant ID"].isin(df_clean["subject_id"]).map({True: "df_clean", False: "not_df_clean"})
        
        if filter2 == 'Patient Count': 
            grouped_df2 = filtered_df2.groupby(["Group", "Material Type"])['Participant ID'].nunique().reset_index(name="Sample Count")
            donut_plot_title = "Donut Plot: Unique Patient by Group and Material Type"
        else:
            grouped_df2 = filtered_df2.groupby(["Group", "Material Type"]).size().reset_index(name="Sample Count")
            donut_plot_title = "Donut Plot: Tubes Count by Group and Material Type"

        # Donut Chart (from grouped_df2)
        donut_df = grouped_df2.groupby("Material Type")["Sample Count"].sum().reset_index()
        donut_df["Percentage"] = donut_df["Sample Count"] / donut_df["Sample Count"].sum()
        donut_df["Percentage Label"] = donut_df["Percentage"].apply(lambda x: f"{x:.1%}")

        # Use a nicer color palette for material types
        color_map = {
            'Plasma': '#4E79A7',      # blue
            'Serum': '#F28E2B',       # orange
            'Biological': '#59A14F',  # green
            'PBMC': '#E15759',        # red
            'Stool': '#B07AA1',       # purple
            'Whole Blood': '#76B7B2'  # teal
        }

        donut = (
            alt.Chart(donut_df)
            .mark_arc(innerRadius=80, outerRadius=180, stroke='white', strokeWidth=2)
            .encode(
            theta=alt.Theta("Sample Count:Q", stack=True),
            color=alt.Color(
                "Material Type:N",
                scale=alt.Scale(domain=list(color_map.keys()), range=list(color_map.values())),
                legend=alt.Legend(title="Material Type", orient="right")
            ),
            tooltip=[
                alt.Tooltip("Material Type:N", title="Material Type"),
                alt.Tooltip("Sample Count:Q", title="Sample Count"),
                alt.Tooltip("Percentage Label:N", title="Percentage")
            ]
            )
            .properties(width=420, height=420, title=donut_plot_title)
        )

    
        
        donut_html = donut.to_html()

        return render_template(
            "material_type.html",
            selected_biorepository=selected_biorepository,
            donut_chart=donut_html,
            grouped_df2=grouped_df2.to_dict(orient='records'),
            donut_plot_title=donut_plot_title,
            cancer_types=cancer_types,
            timepoints=timepoints,
            cancer=cancer,
            timepoint=timepoint,
            filter2=filter2

)

@app.route('/patient/<subject_id>')
def patient_detail(subject_id):
    selected_biorepository = request.form.get('biorepository', '2024')
    biorepo_df = get_biorepo_df(selected_biorepository)

    # Filter all samples for this patient
    patient_samples = biorepo_df[biorepo_df["Participant ID"] == subject_id]

    # Optionally, get patient metadata from df_clean
    all_patient_meta = df_clean[df_clean["subject_id"] == subject_id].to_dict(orient='records')
    patient_meta = {}
    if all_patient_meta:
        for key, value in all_patient_meta[0].items():
            if value != 0 and not pd.isna(value):
                patient_meta[key] = value


    return render_template(
        "patient_detail.html",
        selected_biorepository=selected_biorepository,
        subject_id=subject_id,
        patient_samples=patient_samples.to_dict(orient='records'),
        patient_meta=patient_meta
    )

if __name__ == '__main__':
    app.run(debug=True)
