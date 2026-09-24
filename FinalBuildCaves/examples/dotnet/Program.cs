using System.Runtime.InteropServices;
using System.Text;
using System.Text.Json;

// No MOOCoW reference, database, schema migration or renderer is required.
// On Linux point LD_LIBRARY_PATH at the SDK's installed lib directory.
byte[] request = Encoding.UTF8.GetBytes("{\"world_seed\":42,\"loop_count\":0}");
byte[] error = new byte[1024];
IntPtr result = IntPtr.Zero;
try
{
    int status = Native.Generate(request, (UIntPtr)request.Length, IntPtr.Zero,
        ref result, error, (UIntPtr)error.Length);
    if (status != 0)
        throw new InvalidOperationException($"Native status {status}: {Encoding.UTF8.GetString(error).TrimEnd('\0')}");
    IntPtr bytes = Native.Json(result, out UIntPtr count);
    string json = Marshal.PtrToStringUTF8(bytes, checked((int)count.ToUInt64()))
        ?? throw new InvalidOperationException("Null result");
    using JsonDocument document = JsonDocument.Parse(json);
    Console.WriteLine($"Native C++ generated {document.RootElement.GetProperty("caverns").GetArrayLength()} caverns; schema {document.RootElement.GetProperty("version").GetInt32()}.");
}
finally
{
    Native.Destroy(result);
}

internal static class Native
{
    [DllImport("fbs_caves", EntryPoint = "fbs_caves_generate_json", CallingConvention = CallingConvention.Cdecl)]
    internal static extern int Generate(byte[] request, UIntPtr requestSize, IntPtr options,
        ref IntPtr result, [Out] byte[] error, UIntPtr errorCapacity);
    [DllImport("fbs_caves", EntryPoint = "fbs_caves_result_json", CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr Json(IntPtr result, out UIntPtr size);
    [DllImport("fbs_caves", EntryPoint = "fbs_caves_result_destroy", CallingConvention = CallingConvention.Cdecl)]
    internal static extern void Destroy(IntPtr result);
}
