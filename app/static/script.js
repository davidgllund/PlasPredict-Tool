// DOM Elements
const predictionForm = document.getElementById('predictionForm');
const sequenceSourceRadios = document.querySelectorAll('input[name="sequence_source"]');
const textInputDiv = document.getElementById('textInputDiv');
const fileInputDiv = document.getElementById('fileInputDiv');
const submitBtn = document.getElementById('submitBtn');
const loadingSpinner = document.getElementById('loadingSpinner');
const resultsSection = document.getElementById('resultsSection');
const errorSection = document.getElementById('errorSection');
const newPredictionBtn = document.getElementById('newPredictionBtn');
const closeErrorBtn = document.getElementById('closeErrorBtn');

// Event Listeners
sequenceSourceRadios.forEach(radio => {
    radio.addEventListener('change', updateInputDisplay);
});

// Enforce single-select for isolation_source
document.addEventListener('change', function(e) {
    if (e.target.name === 'isolation_source' && e.target.checked) {
        // Uncheck all other isolation_source checkboxes
        document.querySelectorAll('input[name="isolation_source"]').forEach(checkbox => {
            if (checkbox !== e.target) {
                checkbox.checked = false;
            }
        });
    }
});

predictionForm.addEventListener('submit', handleFormSubmit);
newPredictionBtn.addEventListener('click', resetForm);
closeErrorBtn.addEventListener('click', hideError);

// Initialize
updateInputDisplay();
loadFormOptions();

/**
 * Load form options from API (esp. all inc types from file)
 */
async function loadFormOptions() {
    try {
        const response = await fetch('/api/metadata');
        const data = await response.json();
        
        // Populate inc_types from file
        const incTypeContainer = document.getElementById('inc-type-container');
        if (incTypeContainer && data.inc_types) {
            incTypeContainer.innerHTML = '';
            data.inc_types.forEach(type => {
                const label = document.createElement('label');
                label.innerHTML = `<input type="checkbox" name="inc_type" value="${type}"> ${type}`;
                incTypeContainer.appendChild(label);
            });
        }
    } catch (error) {
        console.error('Error loading form options:', error);
    }
}

/**
 * Update input display based on selected sequence source
 */
function updateInputDisplay() {
    const selectedSource = document.querySelector('input[name="sequence_source"]:checked').value;
    
    if (selectedSource === 'text') {
        textInputDiv.style.display = 'block';
        fileInputDiv.style.display = 'none';
        document.getElementById('sequenceInput').required = true;
        document.getElementById('sequenceFile').required = false;
    } else {
        textInputDiv.style.display = 'none';
        fileInputDiv.style.display = 'block';
        document.getElementById('sequenceInput').required = false;
        document.getElementById('sequenceFile').required = true;
    }
}

/**
 * Handle form submission
 */
async function handleFormSubmit(e) {
    e.preventDefault();
    
    // Hide previous errors
    hideError();
    
    // Validate form
    const validation = validateForm();
    if (!validation.valid) {
        showError(validation.message);
        return;
    }
    
    // Show loading spinner
    loadingSpinner.style.display = 'block';
    submitBtn.disabled = true;
    
    try {
        // Prepare form data
        const formData = new FormData(predictionForm);
        
        // Make API request
        const response = await fetch('/api/predict', {
            method: 'POST',
            body: formData
        });
        
        const result = await response.json();
        
        if (result.success) {
            displayResults(result);
        } else {
            showError(result.error || 'An error occurred during prediction');
        }
    } catch (error) {
        showError(`Network error: ${error.message}`);
    } finally {
        loadingSpinner.style.display = 'none';
        submitBtn.disabled = false;
    }
}

/**
 * Validate form inputs
 */
function validateForm() {
    const sequenceSource = document.querySelector('input[name="sequence_source"]:checked').value;
    const incTypes = document.querySelectorAll('input[name="inc_type"]:checked').length;
    const drugClasses = document.querySelectorAll('input[name="drug_class"]:checked').length;
    const resistanceMechanisms = document.querySelectorAll('input[name="resistance_mechanism"]:checked').length;
    const isolationSources = document.querySelectorAll('input[name="isolation_source"]:checked').length;
    
    if (sequenceSource === 'text') {
        const sequenceText = document.getElementById('sequenceInput').value.trim();
        if (!sequenceText) {
            return { valid: false, message: 'Please provide a FASTA sequence' };
        }
    } else {
        const file = document.getElementById('sequenceFile').files[0];
        if (!file) {
            return { valid: false, message: 'Please select a file' };
        }
    }
    
    // Inc types, drug classes, and resistance mechanisms can be zero or more
    // (no validation required for these)
    
    // Isolation source must be exactly one
    if (isolationSources !== 1) {
        return { valid: false, message: 'Please select exactly one isolation source' };
    }
    
    return { valid: true };
}

/**
 * Display prediction results
 */
function displayResults(result) {
    // Hide form, show results
    predictionForm.style.display = 'none';
    resultsSection.style.display = 'block';
    errorSection.style.display = 'none';
    
    // Populate sequence info
    document.getElementById('resultSeqId').textContent = result.sequence_id;
    document.getElementById('resultSeqLength').textContent = 
        result.sequence_length.toLocaleString();
    
    // Populate conjugation systems if any
    if (result.conjugation_systems && result.conjugation_systems.length > 0) {
        const conjugationCard = document.getElementById('conjugationCard');
        const conjugationSystems = document.getElementById('conjugationSystems');
        conjugationCard.style.display = 'block';
        conjugationSystems.innerHTML = result.conjugation_systems
            .map(system => `<span class="conj-badge">${system}</span>`)
            .join('');
    }
    
    // Populate top prediction
    document.getElementById('topPrediction').textContent = result.top_prediction;
    document.getElementById('topProbability').textContent = 
        result.top_probability.toFixed(3);
    
    // Populate detailed predictions
    displayDetailedPredictions(result.predictions);
    
    // Scroll to results
    resultsSection.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Display detailed prediction bars
 */
function displayDetailedPredictions(predictions) {
    const chart = document.getElementById('predictionsChart');
    chart.innerHTML = '';
    
    // Convert to array and sort by probability (highest first)
    const entries = Object.entries(predictions).sort((a, b) => b[1] - a[1]);
    
    entries.forEach(([label, probability]) => {
        const item = document.createElement('div');
        item.className = 'prediction-item';
        
        const labelDiv = document.createElement('div');
        labelDiv.className = 'prediction-label';
        labelDiv.innerHTML = `
            <span>${label}</span>
            <span>${probability.toFixed(3)}</span>
        `;
        
        const barDiv = document.createElement('div');
        barDiv.className = 'prediction-bar';
        
        const fillDiv = document.createElement('div');
        fillDiv.className = 'prediction-fill';
        fillDiv.style.width = '0%';
        
        barDiv.appendChild(fillDiv);
        
        item.appendChild(labelDiv);
        item.appendChild(barDiv);
        chart.appendChild(item);
        
        // Animate bar fill
        setTimeout(() => {
            fillDiv.style.width = (probability * 100) + '%';
        }, 100);
    });
}

/**
 * Show error message
 */
function showError(message) {
    errorSection.style.display = 'block';
    document.getElementById('errorMessage').textContent = message;
    resultsSection.style.display = 'none';
    predictionForm.style.display = 'block';
    errorSection.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Hide error message
 */
function hideError() {
    errorSection.style.display = 'none';
}

/**
 * Reset form and show input section again
 */
function resetForm() {
    predictionForm.reset();
    resultsSection.style.display = 'none';
    errorSection.style.display = 'none';
    predictionForm.style.display = 'block';
    updateInputDisplay();
    predictionForm.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Load metadata on page load for additional validation/features
 */
async function loadMetadata() {
    try {
        const response = await fetch('/api/metadata');
        const metadata = await response.json();
        // Store in window for potential future use
        window.appMetadata = metadata;
    } catch (error) {
        console.warn('Could not load metadata:', error);
    }
}

// Load metadata on page load
document.addEventListener('DOMContentLoaded', loadMetadata);
