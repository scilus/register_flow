#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["run_bet":"$params.run_bet",
                "quick_registration":"$params.quick_registration",
                "linear_registration":"$params.linear_registration",
                "output_dir":"$params.output_dir",
                "processes_register":"$params.processes_register",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Run Registration to template space"
log.info "============================="
log.info ""

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
log.info "Root: $params.input"
log.info "Template: $params.template"
log.info "Output directory: $params.output_dir"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "Template registration: $params.processes_register"
log.info ""

root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
Channel
    .fromPath("$root/**/*t1.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{in_t1; subjects_for_count}

in_trks = Channel
    .fromFilePairs("$root/**/tractograms/*.trk",
                    size: -1,
                    maxDepth:2) {it.parent.parent.name}

if (!params.template.contains("template_t1.nii.gz")) {
    error "Error ~ Template filename should be named template_t1.nii.gz."}

Channel.fromPath(file(params.template))
    .into{template_for_registration;template_for_transformation_trks;template_for_transformation_metrics; template_check_name}

in_metrics = Channel
    .fromFilePairs("$root/**/metrics/*.nii.gz",
                   size: -1,
                   maxDepth:2) {it.parent.parent.name}

subjects_for_count.count()
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path."}

if (!params.template){
    error "Error ~ Please set the template to use."
}


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
      antsBrainExtraction.sh -d 3 -a ${anat} -e $params.template_t1/t1_template.nii.gz\
            -o bet/ -m $params.template_t1/t1_brain_probability_map.nii.gz -u 0
        scil_image_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz --data_type uint8
        scil_image_math.py multiplication $t1 ${sid}__t1_bet_mask.nii.gz ${sid}__t1_bet.nii.gz
        
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
      """
      scil_apply_transform_to_tractogram.py ${tractogram} ${template} ${transfo} ${tractogram.getSimpleName()}_to_template.trk --inverse
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
      """
      scil_apply_transform_to_tractogram.py ${tractogram} ${template} ${transfo} ${tractogram.getSimpleName()}_to_template.trk --inverse --inverse --in_deformation ${inverse_warp}
      """
}
