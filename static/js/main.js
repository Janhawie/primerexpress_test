/**
 * PCR Primer Design Tool JavaScript
 * Handles UI interactions and client-side functionality
 */

document.addEventListener('DOMContentLoaded', function() {
    // Initialize UI components
    initializeUI();
    
    // Set up event listeners
    setupEventListeners();
    
    // Initialize sequence validator
    if (document.getElementById('sequence_text')) {
        sequenceValidator = new SequenceValidator('sequence_text', 'sequence-validation', 'sequence-stats');
    }
});

/**
 * Initialize UI components
 */
function initializeUI() {
    // Close flash messages when close button is clicked
    document.querySelectorAll('.close-flash').forEach(button => {
        button.addEventListener('click', function() {
            this.closest('.flash').style.display = 'none';
        });
    });
    
    // Modal functionality
    const helpModal = document.getElementById('help-modal');
    if (helpModal) {
        const showHelpBtn = document.getElementById('show-help-btn');
        const closeModal = document.querySelector('.close-modal');
        
        showHelpBtn.addEventListener('click', function(e) {
            e.preventDefault();
            helpModal.style.display = 'flex';
        });
        
        closeModal.addEventListener('click', function() {
            helpModal.style.display = 'none';
        });
        
        // Close modal when clicking outside the content
        window.addEventListener('click', function(e) {
            if (e.target === helpModal) {
                helpModal.style.display = 'none';
            }
        });
    }
    
    // Initialize file upload label
    const fileInput = document.getElementById('sequence_file');
    if (fileInput) {
        fileInput.addEventListener('change', function() {
            const fileName = this.files[0] ? this.files[0].name : 'No file chosen';
            document.getElementById('file-name').textContent = fileName;
        });
    }
}

/**
 * Set up event listeners for the form elements
 */
function setupEventListeners() {
    // Toggle sequence input method
    const sequenceInputMethod = document.getElementById('sequence_input_method');
    if (sequenceInputMethod) {
        sequenceInputMethod.addEventListener('change', toggleSequenceInput);
    }
    
    // Toggle probe parameters
    const includeProbe = document.getElementById('include_probe');
    if (includeProbe) {
        includeProbe.addEventListener('change', toggleProbeParams);
    }
    
    // Use sample sequence button
    const useSampleBtn = document.getElementById('use-sample');
    if (useSampleBtn) {
        useSampleBtn.addEventListener('click', function() {
            document.getElementById('sequence_text').value = sampleSequence;
            sequenceValidator.validate();
        });
    }
    
    // Load preset parameters button
    const loadPresetBtn = document.getElementById('load-preset');
    const presetDropdown = document.querySelector('.preset-dropdown');
    
    if (loadPresetBtn && presetDropdown) {
        loadPresetBtn.addEventListener('click', function() {
            presetDropdown.style.display = presetDropdown.style.display === 'none' ? 'block' : 'none';
        });
        
        // Handle preset selection
        document.querySelectorAll('.preset-option').forEach(option => {
            option.addEventListener('click', function() {
                const presetName = this.dataset.preset;
                if (presets[presetName]) {
                    applyPreset(presets[presetName]);
                    presetDropdown.style.display = 'none';
                }
            });
        });
    }
    
    // Toggle collapsible sections
    document.querySelectorAll('.section-header').forEach(header => {
        header.addEventListener('click', function() {
            toggleCollapsible(this);
        });
    });
    
    // Validate sequence button
    const validateSequenceBtn = document.getElementById('validate-sequence');
    if (validateSequenceBtn) {
        validateSequenceBtn.addEventListener('click', function() {
            sequenceValidator.validate();
        });
    }
    
    // Clear sequence button
    const clearSequenceBtn = document.getElementById('clear-sequence');
    if (clearSequenceBtn) {
        clearSequenceBtn.addEventListener('click', function() {
            document.getElementById('sequence_text').value = '';
            document.getElementById('sequence-validation').textContent = '';
            document.getElementById('sequence-validation').className = 'sequence-feedback';
            document.getElementById('sequence-stats').innerHTML = '';
        });
    }
    
    // Reset parameters button
    const resetParametersBtn = document.getElementById('reset-parameters');
    if (resetParametersBtn) {
        resetParametersBtn.addEventListener('click', function() {
            resetFormParameters();
        });
    }
    
    // Quick design button - uses AJAX
    const quickDesignBtn = document.getElementById('quick-design-btn');
    if (quickDesignBtn) {
        quickDesignBtn.addEventListener('click', function(e) {
            e.preventDefault();
            
            // Validate sequence first
            if (!sequenceValidator.validate()) {
                return false;
            }
            
            // Show loading overlay
            document.getElementById('loading-overlay').style.display = 'flex';
            
            // Get form data
            const formData = getFormData();
            
            // Send AJAX request
            fetch('/api/design', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(formData),
            })
            .then(response => response.json())
            .then(data => {
                // Hide loading overlay
                document.getElementById('loading-overlay').style.display = 'none';
                
                if (data.error) {
                    showQuickResultsError(data.error);
                } else {
                    renderQuickResults(data.results);
                }
            })
            .catch(error => {
                // Hide loading overlay
                document.getElementById('loading-overlay').style.display = 'none';
                showQuickResultsError('An error occurred while processing your request. Please try again.');
                console.error('Error:', error);
            });
        });
    }
    
    // Form submission - show loading overlay
    const primerDesignForm = document.getElementById('primer-design-form');
    if (primerDesignForm) {
        primerDesignForm.addEventListener('submit', function(e) {
            // Show loading overlay
            document.getElementById('loading-overlay').style.display = 'flex';
        });
    }
}

/**
 * Toggle sequence input method (paste or upload)
 */
function toggleSequenceInput() {
    const method = document.getElementById('sequence_input_method').value;
    if (method === 'paste') {
        document.getElementById('paste_input').style.display = 'block';
        document.getElementById('upload_input').style.display = 'none';
    } else {
        document.getElementById('paste_input').style.display = 'none';
        document.getElementById('upload_input').style.display = 'block';
    }
}

/**
 * Toggle probe design parameters visibility
 */
function toggleProbeParams() {
    const includeProbe = document.getElementById('include_probe').checked;
    document.getElementById('probe_parameters').style.display = includeProbe ? 'block' : 'none';
}

/**
 * Toggle collapsible section
 */
function toggleCollapsible(header) {
    const section = header.closest('.collapsible');
    section.classList.toggle('collapsed');
}

/**
 * Apply parameters preset
 */
function applyPreset(preset) {
    // Apply primer parameters
    document.getElementById('primer_opt_size').value = preset.primer_opt_size;
    document.getElementById('primer_min_size').value = preset.primer_min_size;
    document.getElementById('primer_max_size').value = preset.primer_max_size;
    document.getElementById('primer_opt_tm').value = preset.primer_opt_tm;
    document.getElementById('primer_min_tm').value = preset.primer_min_tm;
    document.getElementById('primer_max_tm').value = preset.primer_max_tm;
    document.getElementById('primer_min_gc').value = preset.primer_min_gc;
    document.getElementById('primer_max_gc').value = preset.primer_max_gc;
    document.getElementById('product_min_size').value = preset.product_min_size;
    document.getElementById('product_max_size').value = preset.product_max_size;
    document.getElementById('num_primer_pairs').value = preset.num_primer_pairs;
    
    // Set probe checkbox
    document.getElementById('include_probe').checked = preset.include_probe;
    toggleProbeParams();
    
    // Apply probe parameters if included
    if (preset.include_probe) {
        if (preset.probe_opt_size) document.getElementById('probe_opt_size').value = preset.probe_opt_size;
        if (preset.probe_min_size) document.getElementById('probe_min_size').value = preset.probe_min_size;
        if (preset.probe_max_size) document.getElementById('probe_max_size').value = preset.probe_max_size;
        if (preset.probe_opt_tm) document.getElementById('probe_opt_tm').value = preset.probe_opt_tm;
        if (preset.probe_min_tm) document.getElementById('probe_min_tm').value = preset.probe_min_tm;
        if (preset.probe_max_tm) document.getElementById('probe_max_tm').value = preset.probe_max_tm;
        if (preset.probe_min_gc) document.getElementById('probe_min_gc').value = preset.probe_min_gc;
        if (preset.probe_max_gc) document.getElementById('probe_max_gc').value = preset.probe_max_gc;
    }
}

/**
 * Reset form parameters to defaults
 */
function resetFormParameters() {
    // Reset primer parameters
    document.getElementById('primer_opt_size').value = 20;
    document.getElementById('primer_min_size').value = 18;
    document.getElementById('primer_max_size').value = 27;
    document.getElementById('primer_opt_tm').value = 60.0;
    document.getElementById('primer_min_tm').value = 57.0;
    document.getElementById('primer_max_tm').value = 63.0;
    document.getElementById('primer_min_gc').value = 20.0;
    document.getElementById('primer_max_gc').value = 80.0;
    document.getElementById('product_min_size').value = 100;
    document.getElementById('product_max_size').value = 1000;
    document.getElementById('num_primer_pairs').value = 5;
    
    // Reset probe parameters
    document.getElementById('include_probe').checked = false;
    document.getElementById('probe_parameters').style.display = 'none';
    document.getElementById('probe_opt_size').value = 20;
    document.getElementById('probe_min_size').value = 18;
    document.getElementById('probe_max_size').value = 27;
    document.getElementById('probe_opt_tm').value = 65.0;
    document.getElementById('probe_min_tm').value = 62.0;
    document.getElementById('probe_max_tm').value = 68.0;
    document.getElementById('probe_min_gc').value = 20.0;
    document.getElementById('probe_max_gc').value = 80.0;
}

/**
 * Get form data for AJAX requests
 */
function getFormData() {
    const sequenceText = document.getElementById('sequence_text').value.trim();
    const includeProbe = document.getElementById('include_probe').checked;
    
    const formData = {
        sequence: sequenceText,
        primer_opt_size: parseInt(document.getElementById('primer_opt_size').value),
        primer_min_size: parseInt(document.getElementById('primer_min_size').value),
        primer_max_size: parseInt(document.getElementById('primer_max_size').value),
        primer_opt_tm: parseFloat(document.getElementById('primer_opt_tm').value),
        primer_min_tm: parseFloat(document.getElementById('primer_min_tm').value),
        primer_max_tm: parseFloat(document.getElementById('primer_max_tm').value),
        primer_min_gc: parseFloat(document.getElementById('primer_min_gc').value),
        primer_max_gc: parseFloat(document.getElementById('primer_max_gc').value),
        product_min_size: parseInt(document.getElementById('product_min_size').value),
        product_max_size: parseInt(document.getElementById('product_max_size').value),
        num_primer_pairs: parseInt(document.getElementById('num_primer_pairs').value),
        include_probe: includeProbe
    };
    
    // Add probe parameters if included
    if (includeProbe) {
        formData.probe_opt_size = parseInt(document.getElementById('probe_opt_size').value);
        formData.probe_min_size = parseInt(document.getElementById('probe_min_size').value);
        formData.probe_max_size = parseInt(document.getElementById('probe_max_size').value);
        formData.probe_opt_tm = parseFloat(document.getElementById('probe_opt_tm').value);
        formData.probe_min_tm = parseFloat(document.getElementById('probe_min_tm').value);
        formData.probe_max_tm = parseFloat(document.getElementById('probe_max_tm').value);
        formData.probe_min_gc = parseFloat(document.getElementById('probe_min_gc').value);
        formData.probe_max_gc = parseFloat(document.getElementById('probe_max_gc').value);
    }
    
    return formData;
}

/**
 * Show error message in quick results area
 */
function showQuickResultsError(errorMessage) {
    const quickResults = document.getElementById('quick-results');
    quickResults.style.display = 'block';
    quickResults.innerHTML = `
        <div class="flash flash-error">
            <i class="fas fa-exclamation-circle"></i> ${errorMessage}
        </div>
    `;
}

/**
 * Render quick design results
 */
function renderQuickResults(results) {
    const quickResults = document.getElementById('quick-results');
    quickResults.style.display = 'block';
    
    let html = `
        <h3><i class="fas fa-flask"></i> Quick Design Results</h3>
        <p>Primers designed successfully! Here's a preview of the results:</p>
    `;
    
    // For each sequence
    for (const seqId in results) {
        const seqData = results[seqId];
        
        html += `
            <div class="quick-result-sequence">
                <h4>${seqId}</h4>
                <p>Sequence length: ${seqData.sequence_length} bp</p>
                <p>Primer pairs found: ${seqData.results.length}</p>
                
                <div class="quick-result-primers">
                    <h5>Top Primers:</h5>
                    <table class="primers-table">
                        <thead>
                            <tr>
                                <th>Pair</th>
                                <th>Forward Primer</th>
                                <th>Reverse Primer</th>
                                <th>Size</th>
                            </tr>
                        </thead>
                        <tbody>
        `;
        
        // Show up to 3 primer pairs
        const maxPairs = Math.min(3, seqData.results.length);
        for (let i = 0; i < maxPairs; i++) {
            const pair = seqData.results[i];
            html += `
                <tr>
                    <td>${pair.pair_id}</td>
                    <td class="primer-cell">
                        <div class="primer-sequence">${pair.forward_primer.sequence}</div>
                        <div class="primer-stats">
                            <span>${pair.forward_primer.tm.toFixed(1)}°C</span>
                        </div>
                    </td>
                    <td class="primer-cell">
                        <div class="primer-sequence">${pair.reverse_primer.sequence}</div>
                        <div class="primer-stats">
                            <span>${pair.reverse_primer.tm.toFixed(1)}°C</span>
                        </div>
                    </td>
                    <td>${pair.product_size} bp</td>
                </tr>
            `;
        }
        
        html += `
                        </tbody>
                    </table>
                </div>
            </div>
        `;
    }
    
    html += `
        <div class="quick-result-actions">
            <p>Submit the form to see complete results and download options.</p>
            <button id="close-quick-results" class="btn-secondary">
                <i class="fas fa-times"></i> Close Preview
            </button>
        </div>
    `;
    
    quickResults.innerHTML = html;
    
    // Add event listener to close button
    document.getElementById('close-quick-results').addEventListener('click', function() {
        quickResults.style.display = 'none';
    });
}

/**
 * Sequence Validator Class
 * Validates FASTA sequences and provides statistics
 */
class SequenceValidator {
    constructor(textareaId, validationId, statsId) {
        this.textarea = document.getElementById(textareaId);
        this.validationEl = document.getElementById(validationId);
        this.statsEl = document.getElementById(statsId);
    }
    
    validate() {
        const sequence = this.textarea.value.trim();
        
        // Clear previous validation
        this.validationEl.textContent = '';
        this.validationEl.className = 'sequence-feedback';
        this.statsEl.innerHTML = '';
        
        // Check if empty
        if (!sequence) {
            this.validationEl.textContent = 'Please enter a sequence';
            this.validationEl.className = 'sequence-feedback invalid';
            return false;
        }
        
        // Check if valid FASTA format
        if (!sequence.startsWith('>')) {
            // Try to validate as plain DNA
            const dnaRegex = /^[ATGCNatgcn\s]+$/;
            if (!dnaRegex.test(sequence)) {
                this.validationEl.textContent = 'Invalid sequence format. Please enter a valid DNA sequence or FASTA format';
                this.validationEl.className = 'sequence-feedback invalid';
                return false;
            }
            
            // Valid plain DNA, calculate stats
            const cleanSequence = sequence.replace(/\s/g, '').toUpperCase();
            this.calculateStats(cleanSequence);
            
            this.validationEl.textContent = 'Valid DNA sequence (Plain format)';
            this.validationEl.className = 'sequence-feedback valid';
            return true;
        }
        
        // Parse FASTA format
        const fastaEntries = this.parseFasta(sequence);
        
        if (fastaEntries.length === 0) {
            this.validationEl.textContent = 'Invalid FASTA format';
            this.validationEl.className = 'sequence-feedback invalid';
            return false;
        }
        
        // Calculate stats for each FASTA entry
        for (const entry of fastaEntries) {
            this.calculateStats(entry.sequence, entry.id);
        }
        
        this.validationEl.textContent = `Valid FASTA sequence (${fastaEntries.length} entries)`;
        this.validationEl.className = 'sequence-feedback valid';
        return true;
    }
    
    parseFasta(fastaText) {
        const entries = [];
        const lines = fastaText.split('\n');
        let currentId = '';
        let currentSequence = '';
        
        for (const line of lines) {
            const trimmedLine = line.trim();
            if (!trimmedLine) continue;
            
            if (trimmedLine.startsWith('>')) {
                // If we have a previous entry, save it
                if (currentId && currentSequence) {
                    entries.push({ id: currentId, sequence: currentSequence });
                }
                
                // Start a new entry
                currentId = trimmedLine.substring(1).trim();
                currentSequence = '';
            } else {
                // Add to current sequence, removing any whitespace
                currentSequence += trimmedLine.replace(/\s/g, '').toUpperCase();
            }
        }
        
        // Add the last entry
        if (currentId && currentSequence) {
            entries.push({ id: currentId, sequence: currentSequence });
        }
        
        return entries;
    }
    
    calculateStats(sequence, id = '') {
        const length = sequence.length;
        const a = (sequence.match(/A/g) || []).length;
        const t = (sequence.match(/T/g) || []).length;
        const g = (sequence.match(/G/g) || []).length;
        const c = (sequence.match(/C/g) || []).length;
        const n = (sequence.match(/N/g) || []).length;
        
        const gcContent = ((g + c) / length * 100).toFixed(1);
        
        let statsHtml = '';
        if (id) {
            statsHtml += `<div class="stat-group"><strong>${id}</strong>: `;
        }
        
        statsHtml += `
            <span class="stat">Length: ${length} bp</span>
            <span class="stat">GC: ${gcContent}%</span>
            <span class="stat">A: ${a}</span>
            <span class="stat">T: ${t}</span>
            <span class="stat">G: ${g}</span>
            <span class="stat">C: ${c}</span>
        `;
        
        if (n > 0) {
            statsHtml += `<span class="stat">N: ${n}</span>`;
        }
        
        if (id) {
            statsHtml += '</div>';
        }
        
        this.statsEl.innerHTML += statsHtml;
    }
}