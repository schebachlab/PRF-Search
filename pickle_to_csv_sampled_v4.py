import os
import pickle
import csv
import random
import gc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def sample_background_data(background_category, samples_per_batch):
    """
    Samples data from pickle files for a given background category and saves to a CSV file.

    Parameters:
        background_category (str): The category of background distances ('All', 'Canonical', 'Alternative').
        samples_per_batch (list): List of sample sizes to retrieve from each batch file.
    """
    # List all batch files for the category
    batch_files = [f for f in os.listdir() if f.startswith(f"background_distances_{background_category}_batch_") and f.endswith('.pkl')]

    if not batch_files:
        print(f"No background files found for category {background_category}.")
        return []

    # Sort the batch files to process them in order
    batch_files.sort()

    # Check if the number of batches matches the number of samples_per_batch entries
    if len(batch_files) != len(samples_per_batch):
        print(f"Number of batches ({len(batch_files)}) does not match number of sample sizes ({len(samples_per_batch)}) for category {background_category}.")
        return []

    # Initialize list to hold sampled data
    sampled_data = []

    # Process each batch file
    for idx, batch_file in enumerate(batch_files):
        print(f"Processing {batch_file}")
        with open(batch_file, 'rb') as f:
            # Load data
            data = pickle.load(f)
            batch_size = len(data)
            samples_in_this_batch = samples_per_batch[idx]
            # Sample data
            if batch_size > samples_in_this_batch:
                batch_sample = random.sample(data, samples_in_this_batch)
            else:
                batch_sample = data  # If data is less than samples_in_this_batch, use all data
            # Add sampled data to the list
            sampled_data.extend(batch_sample)
            # Clean up
            del data, batch_sample
            gc.collect()

    # Output CSV file name (optional)
    csv_file_name = f"sampled_v4_background_distances_{background_category}.csv"

    # Write sampled data to CSV (optional)
    with open(csv_file_name, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(['Distance_from_TMD', 'Category'])
        for distance in sampled_data:
            csv_writer.writerow([distance, background_category])

    print(f"Sampled CSV file {csv_file_name} created successfully.")

    return sampled_data

def load_actual_data(actual_csv_file, category):
    """
    Loads actual distances for a given category from the CSV file without sampling.

    Parameters:
        actual_csv_file (str): Path to the CSV file containing actual distances.
        category (str): The category to load (e.g., 'strict_all').

    Returns:
        DataFrame: A pandas DataFrame containing the actual distances for the specified category.
    """
    # Read the actual distances CSV file
    df = pd.read_csv(actual_csv_file)

    # Filter for the specified category
    df_category = df[df['Category'] == category]

    if df_category.empty:
        print(f"No data found for category '{category}' in the actual distances CSV file.")
        return None

    return df_category

def plot_kde(actual_df, background_data, actual_category, background_category, output_file, title, xlim=None):
    """
    Plots a KDE of the actual and background distances.

    Parameters:
        actual_df (DataFrame): DataFrame containing actual distances.
        background_data (list): List of sampled background distances.
        actual_category (str): Category name for the actual data.
        background_category (str): Category name for the background data.
        output_file (str): Path to save the output plot image.
        title (str): Title of the plot.
        xlim (tuple): Optional x-axis limits as a tuple (xmin, xmax).
    """
    # Create DataFrames for actual and background data
    actual_data = pd.DataFrame({
        'Distance_from_TMD': actual_df['Distance_from_TMD'],
        'Category': actual_category
    })

    background_df = pd.DataFrame({
        'Distance_from_TMD': background_data,
        'Category': background_category
    })

    # Initialize the plot
    plt.figure(figsize=(10, 6))

    # Plot KDE for background data
    sns.kdeplot(
        data=background_df,
        x='Distance_from_TMD',
        label=f'Background ({background_category})',
        color='blue',
        linestyle='--'
    )

    # Plot KDE for actual data
    sns.kdeplot(
        data=actual_data,
        x='Distance_from_TMD',
        label=f'Actual ({actual_category})',
        color='red'
    )

    # Set plot title and labels
    plt.title(title)
    plt.xlabel('Distance from Nearest Upstream TMD (codons)')
    plt.ylabel('Density')

    # Set x-axis limits if provided
    if xlim is not None:
        plt.xlim(xlim)

    # Show legend
    plt.legend(title='Category')

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"KDE plot saved to {output_file}")

# Main script

# Mapping of actual categories to background categories
category_mapping = {
    'strict_all': 'All',
    'loose_all': 'All',
    'strict_canonical': 'Canonical',
    'loose_canonical': 'Canonical',
    'strict_alternative': 'Alternative',
    'loose_alternative': 'Alternative'
}

# Load actual data for all categories without sampling
actual_csv_file = 'motif_search_output_v57_production_version_final_v55_slipsite_tmd_distances_strict_and_loose_slipsite_distribution.csv'
actual_categories = list(category_mapping.keys())
actual_data = {}  # Dictionary to store actual data DataFrames
actual_data_sizes = {}  # Dictionary to store actual data sizes per actual category

for actual_category in actual_categories:
    actual_df = load_actual_data(actual_csv_file, actual_category)
    if actual_df is not None:
        actual_data[actual_category] = actual_df
        actual_data_sizes[actual_category] = len(actual_df)
    else:
        print(f"No actual data available for category '{actual_category}'.")
        actual_data_sizes[actual_category] = 0

# Compute desired sample sizes per background category
background_categories = set(category_mapping.values())
desired_sample_sizes = {}

for background_category in background_categories:
    # Get the actual categories mapping to this background category
    actual_categories_for_background = [act_cat for act_cat, back_cat in category_mapping.items() if back_cat == background_category]
    # Sum the sizes of actual data for these categories
    total_actual_size = sum(actual_data_sizes[act_cat] for act_cat in actual_categories_for_background)
    desired_sample_sizes[background_category] = total_actual_size

# Compute samples per batch for each background category
samples_per_batch_per_category = {}
batches_per_category = 8  # Since there are 8 batches per category

for background_category in background_categories:
    sample_size = desired_sample_sizes[background_category]
    base_samples_per_batch = sample_size // batches_per_category
    remainder = sample_size % batches_per_category

    # Initialize list with base_samples_per_batch for each batch
    samples_per_batch = [base_samples_per_batch] * batches_per_category

    # Distribute the remainder
    for i in range(remainder):
        samples_per_batch[i] += 1

    samples_per_batch_per_category[background_category] = samples_per_batch

# Sample background data for each background category using samples per batch
background_samples = {}  # Dictionary to store sampled background data
for background_category in background_categories:
    samples_per_batch = samples_per_batch_per_category[background_category]
    sampled_data = sample_background_data(background_category, samples_per_batch=samples_per_batch)
    background_samples[background_category] = sampled_data

# Plot KDE for each actual category mapped to its background category
for actual_category, background_category in category_mapping.items():
    actual_df = actual_data.get(actual_category)
    background_data = background_samples.get(background_category)

    if actual_df is not None and background_data is not None:
        # Full-range plot
        output_plot_full = f'distance_distribution_kde_{actual_category}_full_range.png'
        title_full = f'Distribution of Distances from TMD ({actual_category}) - Full Range'
        plot_kde(
            actual_df,
            background_data,
            actual_category,
            background_category,
            output_plot_full,
            title_full
        )

        # Zoomed-in plot over range -12 to 12
        output_plot_zoom = f'distance_distribution_kde_{actual_category}_zoomed.png'
        title_zoom = f'Distribution of Distances from TMD ({actual_category}) - Zoomed In (-12 to 12)'
        plot_kde(
            actual_df,
            background_data,
            actual_category,
            background_category,
            output_plot_zoom,
            title_zoom,
            xlim=(-12, 12)  # Set x-axis limits for zoomed-in plot
        )
    else:
        print(f"Data not available for plotting category '{actual_category}'.")
