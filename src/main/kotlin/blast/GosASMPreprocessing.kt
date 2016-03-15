package blast

import com.beust.klaxon.*
import com.google.common.collect.Lists
import java.io.*
import java.nio.file.Files
import java.util.*

val FILE_PATH = "/Users/aganezov/Desktop/Assembly/Allpaths-LG/"

fun main(args: Array<String>) {

}

private fun analyze_contigs_alignment() {
    val contigs = read_fasta(file_path = FILE_PATH, file_name = "genome.ctg.fasta").map {
        it.first to it.second.length
    }
    val chr_14 = read_fasta(file_path = FILE_PATH, file_name = "chr14_human_seq.fa").map {
        it.first to it.second.length
    }.first()

    println("Overall length of chr14:   ${chr_14.second}")
    println("Number of initial contigs: ${contigs.size}")
    println("Overall length of contigs: ${contigs.map { it.second }.sum()}")


    val contigsAlignment = read_mapped_contigs(file_path = FILE_PATH, file_name = "genome.ctg.positions.txt")
    val alignedContigsWithMultipleEntries = contigsAlignment.map { it.fragmentName }.groupBy { it }.filter { it.value.size > 1 }


    val contigsMap = contigs.map { it.first to it.second }.toMap()

    val mappedContigs = contigsAlignment
            .filter { it.fragmentName !in alignedContigsWithMultipleEntries.keys }
            .map {
                FragmentAlignment(
                        fragmentName = it.fragmentName,
                        fragmentStart = it.fragmentStart,
                        fragmentEnd = it.fragmentEnd,
                        fragmentStrand = it.fragmentStrand,
                        fragmentLength = contigsMap[it.fragmentName]?.toLong() ?: -1,

                        subjectName = it.subjectName,
                        subjectStart = it.subjectStart,
                        subjectEnd = it.subjectEnd,
                        subjectStrand = it.subjectStrand,

                        alignmentLength = it.alignmentLength,
                        gapsCnt = it.gapsCnt,
                        identity = it.identity
                )
            }
    println("Number of non duplicated contigs: ${mappedContigs.size}")
    println("Overall coverage by non duplicated contigs ${mappedContigs.map { it.fragmentLength }.sum()}")
    val not_aligned_contigs = mappedContigs.filter { it.fragmentStart == -1L }
    println("Not aligned contigs cnt ${not_aligned_contigs.size}. Their overall coverage ${not_aligned_contigs.map { it.fragmentLength }.sum() }")
    val fullyAlignedContigs = mappedContigs.filter { it.fragmentStart == 1L && it.fragmentEnd in arrayOf(it.fragmentLength, it.fragmentLength - 1, it.fragmentLength + 1) }
    println("Fully aligned contigs cnt ${fullyAlignedContigs.size}. Their overall coverage ${fullyAlignedContigs.map { it.fragmentLength }.sum()}")
    val partiallyAlignedContigs = mappedContigs.filter { it.fragmentStart != 1L || it.fragmentEnd !in arrayOf(it.fragmentLength, it.fragmentLength - 1, it.fragmentLength + 1) }
    val averageMappedPortion = partiallyAlignedContigs.map { Math.abs(it.fragmentEnd - it.fragmentStart) * 100.0 / it.fragmentLength }.average()
    println("Partially aligned contigs cnt ${partiallyAlignedContigs.size}\n" +
            "Their mapped overall coverage ${partiallyAlignedContigs.map { Math.abs(it.fragmentEnd - it.fragmentStart) }.sum()}\n" +
            "Their non-mapped coverage ${partiallyAlignedContigs.map { it.fragmentLength - Math.abs(it.fragmentEnd - it.fragmentStart) }.sum()}\n" +
            "Average mapped portion $averageMappedPortion")
    val sortedMappedContig = mappedContigs.sortedBy { Math.min(it.subjectStart, it.subjectEnd) }
    val sortedFullyMappedContig = fullyAlignedContigs.sortedBy { Math.min(it.subjectStart, it.subjectEnd) }
    val fullyAlignedOverlapping = sortedFullyMappedContig.dropLast(1).zip(sortedFullyMappedContig.drop(1)).filter { Math.max(it.first.subjectStart, it.first.subjectEnd) > Math.min(it.second.subjectStart, it.second.subjectEnd) }
    val mappedOverlapping = sortedMappedContig.dropLast(1).zip(sortedMappedContig.drop(1)).filter { Math.max(it.first.subjectStart, it.first.subjectEnd) > Math.min(it.second.subjectStart, it.second.subjectEnd) }
    println("Number of mapped overlapping contigs : ${mappedOverlapping.size}. ${mappedOverlapping.map { it.first.fragmentName to it.second.fragmentName }}")
    println("Number of fully aligned overlapping contigs: ${fullyAlignedOverlapping.size}. ${fullyAlignedOverlapping.map { it.first.fragmentName to it.second.fragmentName }}")
}

private fun obtain_contig_alignments() {
    var flag: Boolean
    val alignmentResultsWriter = BufferedWriter(FileWriter("/Users/aganezov/Desktop/Assembly/Allpaths-LG/genome.ctg.positions.txt"))
    alignmentResultsWriter.write("" +
            "fragment_name\t" +
            "fragment_start\t" +
            "fragment_end\t" +
            "fragment_strand\t" +

            "subject_name\t" +
            "subject_start\t" +
            "subject_end\t" +
            "subject_strand\t" +

            "alignment_length\t" +
            "identity\t" +
            "gaps_cnt")
    alignmentResultsWriter.newLine()
    alignmentResultsWriter.flush()
    val alreadyProcessedFiles = mutableListOf<String>();

    try {
        do {
            try {
                flag = run_blast_for_fasta_files(
                        directory_with_query_files = "/Users/aganezov/Desktop/Assembly/Allpaths-LG/splitted_fasta",
                        subject_file_full_path = "/Users/aganezov/Desktop/Assembly/Allpaths-LG/chr14_human_seq.fa",
                        result_destination_directory = "/Users/aganezov/Desktop/Assembly/Allpaths-LG/alignments",
                        alreadyProcessedFiles = alreadyProcessedFiles,
                        resultsWriter = alignmentResultsWriter
                )
            } catch(e: IllegalThreadStateException) {
                println(e.message)
                println(e.stackTrace.toString())
                flag = false
            }
        } while (!flag)
    } finally {
        alignmentResultsWriter.flush()
        alignmentResultsWriter.close()
    }
}

fun read_mapped_contigs(file_path: String, file_name: String): List<FragmentAlignment> {
    val reader = BufferedReader(FileReader(file_path + File.separator + file_name))
    reader.readLine()
    return reader.readLines().map {
        val data = it.trim().split("\t")
        FragmentAlignment(
                fragmentName = data[0],
                fragmentStart = data[1].toLong(),
                fragmentEnd = data[2].toLong(),
                fragmentStrand = data[3].toInt(),
                subjectName = data[4],
                subjectStart = data[5].toLong(),
                subjectEnd = data[6].toLong(),
                subjectStrand = data[7].toInt(),
                alignmentLength = data[8].toLong(),
                identity = data[9].toLong(),
                gapsCnt = data[10].toLong()
        )
    }
}

fun parseJson(name: String): Any {
    return Parser().parse(FileInputStream(name))!!
}

private fun convertStrand(strandString: String): Int {
    return when (strandString) {
        "Minus" -> -1
        "Plus" -> 1
        else -> 0
    }
}

private fun get_real_contig_alignment(jsonBlastResultFilePath: String): List<FragmentAlignment> {
    val topLevel = parseJson(jsonBlastResultFilePath) as JsonObject
    val alignmentObject = topLevel.obj("BlastOutput2")!!.obj("report")!!.obj("results")!!.array<JsonObject>("bl2seq")!!.first()
    val queryLength = alignmentObject.long("query_len")
    return alignmentObject.array<JsonObject>("hits")?.firstOrNull()?.array<JsonObject>("hsps")
            //            .filter {
            //                it.long("query_from") as Long == 1L && it.long("query_to") as Long == queryLength
            //            }
            ?.sortedByDescending {
                it.long("identity")
            }
            ?.map {
                FragmentAlignment(
                        fragmentName = alignmentObject.string("query_title") ?: "",
                        fragmentLength = queryLength?.toLong() ?: -1L,
                        subjectName = "chr14",

                        fragmentStrand = convertStrand(it.string("query_strand") ?: ""),
                        subjectStrand = convertStrand(it.string("hit_strand") ?: ""),

                        fragmentStart = it.long("query_from") ?: -1L,
                        fragmentEnd = it.long("query_to") ?: -1L,

                        subjectStart = it.long("hit_from") ?: -1L,
                        subjectEnd = it.long("hit_to") ?: -1L,

                        gapsCnt = it.long("gaps") ?: -1L,
                        identity = it.long("identity") ?: -1L,
                        alignmentLength = it.long("align_len") ?: -1
                )
            }?.toList() ?: listOf(FragmentAlignment(fragmentName = alignmentObject.string("query_title") ?: "----"))
}

private fun run_blast_for_fasta_files(directory_with_query_files: String,
                                      subject_file_full_path: String,
                                      result_destination_directory: String,
                                      alreadyProcessedFiles: MutableList<String>,
                                      resultsWriter: BufferedWriter): Boolean {
    val runtime = Runtime.getRuntime()
    File(directory_with_query_files).listFiles()
            .filter { it.isFile && it.name.endsWith("fasta") }
            .filter { it.name !in alreadyProcessedFiles }
            .map {
                val destination_file_full_path = result_destination_directory + File.separator + it.name.substringBefore(".")

                val command = "/usr/local/ncbi/blast/bin/blastn -query ${it.absolutePath} -subject $subject_file_full_path " +
                        "-out $destination_file_full_path " +
                        "-outfmt 13 " +
                        "-evalue 0.001"
                val p = runtime.exec(command)
                val error_messages = BufferedReader(InputStreamReader(p.errorStream)).readLines()
                val output_message = BufferedReader(InputStreamReader(p.inputStream)).readLines()
                if (error_messages.size > 0 || output_message.size > 0) {
                    println("***************************************************************************")
                    println("***************************************************************************")
                    println("During execution on ${it.name}, following blast output was produced")
                    println("Errors:")

                    println(error_messages.joinToString("\n"))
                    println("Output:")
                    println(output_message.joinToString("\n"))
                }

                if (p.exitValue() == 0) {
                    alreadyProcessedFiles.add(it.name)
                    println("Successfully finished executing blastn on ${it.name}")
                    File(result_destination_directory).listFiles()
                            .filter { it.name.endsWith(".json") }
                            .map {
                                val goodAlignmentResults = get_real_contig_alignment(jsonBlastResultFilePath = it.absolutePath)
                                        .groupBy { it.identity }
                                        .maxBy { it.key }?.value ?: mutableListOf(FragmentAlignment())
                                if (goodAlignmentResults.first().fragmentName != "")
                                    println("Successfully obtained results for ${goodAlignmentResults.first().fragmentName}. " +
                                            "Number of perfect alignments: ${goodAlignmentResults.size}")
                                else
                                    println("Failed to obtain results for ${it.name}. Writing an empty alignment")
                                goodAlignmentResults.forEach {
                                    resultsWriter.write("" +
                                            "${it.fragmentName}\t" +
                                            "${it.fragmentStart}\t" +
                                            "${it.fragmentEnd}\t" +
                                            "${it.fragmentStrand}\t" +

                                            "${it.subjectName}\t" +
                                            "${it.subjectStart}\t" +
                                            "${it.subjectEnd}\t" +
                                            "${it.subjectStrand}\t" +

                                            "${it.alignmentLength}\t" +
                                            "${it.identity}\t" +
                                            "${it.gapsCnt}")
                                    resultsWriter.newLine()
                                    resultsWriter.flush()
                                }
                            }
                    File(result_destination_directory).listFiles().forEach {
                        Files.delete(it.toPath())
                    }
                    if (File(result_destination_directory).listFiles().size > 0) {
                        System.exit(-1)
                    }
                    println("Successfully retrieved all alignment results")
                    println()
                } else {
                    println("Failed to run on ${it.name}")
                    println()
                }
            }
    return true;


}


private fun read_fasta(file_path: String, file_name: String): MutableList<Pair<String, String>> {
    val result = ArrayList<Pair<String, String>>()
    val reader = BufferedReader(FileReader(File(file_path + File.separator + file_name)))
    var current_fragment_name = ""
    var current_fragment_sequence = StringBuilder()
    var line = reader.readLine()
    while (line != null) {
        line = line.trim()
        if (line.startsWith(">")) {
            if (current_fragment_sequence.length > 0 && current_fragment_name.length > 0) {
                result.add(current_fragment_name to current_fragment_sequence.toString())
                current_fragment_sequence = StringBuilder()
            }
            current_fragment_name = line.substring(1)
        } else {
            current_fragment_sequence.append(line)
        }
        line = reader.readLine()
    }
    if (current_fragment_sequence.length > 0) {
        result.add(current_fragment_name to current_fragment_sequence.toString())
    }
    return result
}

private fun split_contigs_into_chunks(file_path: String, file_name: String, destination_path: String,
                                      destination_prefix: String) {
    val data = read_fasta(file_path = file_path, file_name = file_name)
    val val_split_set = Lists.partition(data, 10)
    val_split_set.mapIndexed { index, data ->
        val writer = BufferedWriter(FileWriter(destination_path + File.separator +
                destination_prefix + index.toString() + ".fasta"))
        data.map {
            writer.write(">" + it.first.toString())
            writer.newLine()
            writer.write(it.second)
            writer.newLine()
        }
        writer.flush()
        writer.close()
    }
}

private fun extract_14_chromosome_fasta(file_path: String, file_name: String, destination_file_name: String) {
    val reader = BufferedReader(FileReader(File(file_path + File.separator + file_name)))
    val writer = BufferedWriter(FileWriter(File(file_path + File.separator + destination_file_name)))
    var line: String? = reader.readLine()
    var chr_14_flag = false
    while (line != null) {
        if (line.trim().startsWith(">")) {
            println("processing chr ${line.trim()}")
            chr_14_flag = line.trim().contains("chr14")
            if (line.trim().contains("chr15")) break
        }
        if (chr_14_flag) {
            writer.write(line)
        }
        line = reader.readLine()
    }
    writer.flush()
}