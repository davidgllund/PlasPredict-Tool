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
 * Load form options from API
 */
async function loadFormOptions() {
    try {
        const response = await fetch('/api/metadata');
        const data = await response.json();
        // Metadata loaded for future use
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
    
    // Populate detected features
    if (result.detected_features) {
        displayDetectedFeatures(result.detected_features);
    }
    
    // Populate top prediction
    document.getElementById('topPrediction').textContent = result.top_prediction;
    document.getElementById('topProbability').textContent = 
        result.top_probability.toFixed(3);
    
    // Populate detailed predictions (with fixed order)
    displayDetailedPredictions(result.predictions);
    
    // Scroll to results
    resultsSection.scrollIntoView({ behavior: 'smooth' });
}

/**
 * Display detected features from sequence analysis
 */
function displayDetectedFeatures(features) {
    const detectedFeaturesCard = document.getElementById('detectedFeaturesCard');
    
    if (!detectedFeaturesCard) return;
    
    let hasFeatures = false;
    let featureHtml = '';
    
    // Conjugation systems
    if (features.conjugation_systems && features.conjugation_systems.length > 0) {
        hasFeatures = true;
        featureHtml += '<div class="feature-category"><strong>Conjugation Systems:</strong> ';
        featureHtml += features.conjugation_systems
            .map(system => `<span class="feature-badge">${system}</span>`)
            .join(' ');
        featureHtml += '</div>';
    }
    
    // Inc types
    if (features.inc_types && features.inc_types.length > 0) {
        hasFeatures = true;
        featureHtml += '<div class="feature-category"><strong>Incompatibility Types:</strong> ';
        featureHtml += features.inc_types
            .map(type => `<span class="feature-badge">${type}</span>`)
            .join(' ');
        featureHtml += '</div>';
    }
    
    // Drug classes
    if (features.drug_classes && features.drug_classes.length > 0) {
        hasFeatures = true;
        featureHtml += '<div class="feature-category"><strong>Antibiotic Classes:</strong> ';
        featureHtml += features.drug_classes
            .map(drug => `<span class="feature-badge">${drug}</span>`)
            .join(' ');
        featureHtml += '</div>';
    }
    
    // Resistance mechanisms
    if (features.resistance_mechanisms && features.resistance_mechanisms.length > 0) {
        hasFeatures = true;
        featureHtml += '<div class="feature-category"><strong>Resistance Mechanisms:</strong> ';
        featureHtml += features.resistance_mechanisms
            .map(mech => `<span class="feature-badge">${mech}</span>`)
            .join(' ');
        featureHtml += '</div>';
    }
    
    if (hasFeatures) {
        detectedFeaturesCard.style.display = 'block';
        document.getElementById('detectedFeatures').innerHTML = featureHtml;
    } else {
        detectedFeaturesCard.style.display = 'none';
    }
}

/**
 * Display detailed prediction bars with fixed order
 */
function displayDetailedPredictions(predictions) {
    const chart = document.getElementById('predictionsChart');
    chart.innerHTML = '';
    
    // Define fixed order (top to bottom as requested)
    const fixedOrder = [
        '≥ Phylum',
        'Order/Class',
        'Genus/Family',
        '≤ Species'
    ];
    
    // Display in fixed order
    fixedOrder.forEach(label => {
        const probability = predictions[label];
        if (probability === undefined) return;
        
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
