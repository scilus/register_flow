#!/usr/bin/env nextflow

params.input = false
params.template = false

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["input":"$params.input",
                "template":"$params.template",
                "run_bet":"$params.run_bet",
                "quick_registration":"$params.quick_registration",
                "linear_registration":"$params.linear_registration",
                "trk_keep_invalid":"$params.trk_keep_invalid",
                "trk_cut_invalid":"$params.trk_cut_invalid",
                "trk_remove_invalid":"$params.trk_remove_invalid",
                "output_dir":"$params.output_dir",
                "processes_register":"$params.processes_register",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

root = file(params.input)

if (!params.input){
  error "Error ~ Please set the input directory with --input."}
else {
  Channel
      .fromPath("$root/**/*t1.nii.gz",
                      maxDepth:1, checkIfExists: true)
      .map{[it.parent.name, it]}
      .into{in_t1; subjects_for_count}

  in_metrics = Channel
    .fromFilePairs("$root/**/metrics/*.nii.gz",
                   size: -1,
                   maxDepth:2) {it.parent.parent.name}

  in_trks = Channel
    .fromFilePairs("$root/**/tractograms/*.trk",
                    size: -1,
                    maxDepth:2) {it.parent.parent.name}
}

if (!params.trk_cut_invalid && !params.trk_keep_invalid && !params.trk_remove_invalid){
  log.warn "No option is set to handle invalid streamlines after registration. Default is to remove them."
  trk_remove_invalid = true
}
else{
  trk_remove_invalid = false
}

if ((params.trk_keep_invalid && params.trk_cut_invalid) ||
    (params.trk_keep_invalid && params.trk_remove_invalid) ||
    (params.trk_cut_invalid && params.trk_remove_invalid)){
  log.error "Only one option is allowed to handle invalid streamlines after registration."
  error "Please set only one of the following options: --trk_keep_invalid, --trk_cut_invalid, --trk_remove_invalid."
}

subjects_for_count.count()
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path."}


if (!params.template){
  error "Error ~ Please set the template with --template."}
else
{
  Channel.fromPath(file(params.template), checkIfExists: true)
    .into{template_for_registration;template_for_transformation_trks;template_for_transformation_metrics; template_check_name}
}

log.info ""
log.info "Run Registration to template space"
log.info "=================================="
log.info ""

log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "[Inputs]"
log.info "Input folder: $params.input"
log.info "Template: $params.template"
log.info "Output directory: $params.output_dir"
log.info ""

log.info "[Parameters]"
log.info "Run BET: $params.run_bet"
log.info "Quick registration: $params.quick_registration"
log.info "Linear registration: $params.linear_registration"
log.info ""

log.info "[TRK Parameters]"
log.info "Keep invalid: $params.trk_keep_invalid"
log.info "Cut invalid: $params.trk_cut_invalid"
log.info "Remove invalid: $trk_remove_invalid"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "Template registration: $params.processes_register"
log.info ""

in_t1.combine(template_for_registration).set{anats_for_registration}

process Register_T1_to_template {
    cpus params.processes_register
    memory '2 GB'
    publishDir = params.registration

    input:
    set sid, file(anat), file(template) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1Warp.nii.gz", "${sid}__output1InverseWarp.nii.gz" into nonlinear_transformation_for_trks, nonlinear_transformation_for_metrics optional true
    set sid, "${sid}__output0GenericAffine.mat" into linear_transformation_for_trks, linear_transformation_for_metrics optional true
    set sid, "${sid}__t1_template_space.nii.gz" into t1_to_template
    file "${sid}__t1_bet_mask.nii.gz" optional true
    file "${sid}__t1_bet.nii.gz" optional true

    script:
    if (params.run_bet){
    """
      antsBrainExtraction.sh -d 3 -a ${anat} -e $params.template_t1_for_bet/t1_template.nii.gz\
            -o bet/ -m $params.template_t1_for_bet/t1_brain_probability_map.nii.gz -u 0
        scil_volume_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz --data_type uint8
        scil_volume_math.py multiplication $t1 ${sid}__t1_bet_mask.nii.gz ${sid}__t1_bet.nii.gz
        
      ${params.script_registration} -d 3 -m ${sid}__t1_bet.nii.gz -f ${template} -n ${task.cpus} -o "${sid}__output" -t ${params.transfo}
      mv ${sid}__outputWarped.nii.gz ${sid}__t1_template_space.nii.gz
    """
    }
    else{
    """
      ${params.script_registration} -d 3 -m ${anat} -f ${template} -n ${task.cpus} -o "${sid}__output" -t ${params.transfo}
      mv ${sid}__outputWarped.nii.gz ${sid}__t1_template_space.nii.gz
    """
    }
}

Channel.empty()
  .into{nonlinear_metrics_transformation_for_metrics; linear_metrics_transformation_for_metrics}

if(params.linear_registration){
  in_metrics
      .transpose()
      .combine(template_for_transformation_metrics)
      .combine(linear_transformation_for_metrics, by: 0)
      .set{linear_metrics_transformation_for_metrics}
}
else{
  in_metrics
    .transpose()
    .combine(template_for_transformation_metrics)
    .combine(nonlinear_transformation_for_metrics, by: 0)
    .set{nonlinear_metrics_transformation_for_metrics}
}

process Linear_Registration_Metrics_to_template {
    cpus 1
    publishDir = params.registration_metrics

    input:
    set sid, file(metric), file(template), file(transfo) from linear_metrics_transformation_for_metrics

    output:
    file "*_to_template.nii.gz"

    script:
    if (params.linear_registration)
    """
    if [[ "$metric" == *"nufo"* ]]; then
      antsApplyTransforms -d 3 -i $metric -r $template -t $transfo -o ${metric.getSimpleName()}_to_template.nii.gz -n NearestNeighbor
    else
      antsApplyTransforms -d 3 -i $metric -r $template -t $transfo -o ${metric.getSimpleName()}_to_template.nii.gz
    fi
    """
}

process NonLinear_Registration_Metrics_to_template {
    cpus 1
    publishDir = params.registration_metrics

    input:
    set sid, file(metric), file(template), file(transfo), file(warp), file(inverse_warp) from nonlinear_metrics_transformation_for_metrics

    output:
    file "*_to_template.nii.gz"

    script:
      """
      antsApplyTransforms -d 3 -i $metric -r $template -t $warp $transfo -o ${metric.getSimpleName()}_to_template.nii.gz
      """
}

Channel.empty()
  .into{linear_trks_transformation_for_trks; nonlinear_trks_transformation_for_trks}

if(params.linear_registration){
in_trks
    .transpose()
    .combine(template_for_transformation_trks)
    .combine(linear_transformation_for_trks, by: 0)
    .set{linear_trks_transformation_for_trks}
}
else{
in_trks
    .transpose()
    .combine(template_for_transformation_trks)
    .combine(nonlinear_transformation_for_trks, by: 0)
    .set{nonlinear_trks_transformation_for_trks}
}

process Linear_Registration_Tractograms_to_template {
    cpus 1
    publishDir = params.registration_trks

    input:
    set sid, file(tractogram), file(template), file(transfo) from linear_trks_transformation_for_trks

    output:
    file "*_to_template.trk"

    script:
      def option = trk_remove_invalid ? "--remove_invalid" : ""
      if (params.trk_cut_invalid) {
        option = "--cut_invalid"
      }
      else if (params.trk_keep_invalid) {
        option = "--keep_invalid"
      }
      """
      scil_tractogram_apply_transform.py $option ${tractogram} ${template} ${transfo} ${tractogram.getSimpleName()}_to_template.trk --inverse
      """
}

process NonLinear_Registration_Tractograms_to_template {
    cpus 1
    publishDir = params.registration_trks

    input:
    set sid, file(tractogram), file(template), file(transfo), file(warp), file(inverse_warp) from nonlinear_trks_transformation_for_trks

    output:
    file "*_to_template.trk"

    script:
      def option = trk_remove_invalid ? "--remove_invalid" : ""
      if (params.trk_cut_invalid) {
        option = "--cut_invalid"
      }
      else if (params.trk_keep_invalid) {
        option = "--keep_invalid"
      }
      """
      scil_tractogram_apply_transform.py $option ${tractogram} ${template} ${transfo} ${tractogram.getSimpleName()}_to_template.trk --inverse --inverse --in_deformation ${inverse_warp}
      """
}
